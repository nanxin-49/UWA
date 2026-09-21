function outputs = generate_bellhop_100m_source1m_visuals_vertical(overrides)
% A/E/R illustration of a physical source 1 m above a 100 m seabed.
% Isolated reporting entry: no PE, source-limit audit, or retired data.
if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
cfg=struct('bellhop_exe',resolve_bellhop_exe_vertical(), ...
    'output_dir',fullfile(root,'results','validation','bellhop_100m_rough_surface','source_1m'), ...
    'water_depth_m',100,'source_clearance_m',1,'rx_depth_m',3, ...
    'geometry_angle_deg',89.5,'frequency_hz',4000,'datum_shift_m',1, ...
    'beam_count',10001,'fan_beams',1001,'angle_limits_deg',[-89.9999 -60], ...
    'range_box_m',12,'fan_range_m',120,'step_m',0.02,'surface_dx_m',0.01, ...
    'c0_mps',1500,'bottom_c_mps',1800,'bottom_rho_g_cm3',2, ...
    'bottom_alpha_db_lambda',0.5,'show_figures',false,'render_only',false);
names=fieldnames(overrides);
for ii=1:numel(names)
    assert(isfield(cfg,names{ii}),'Unknown override: %s',names{ii});
    cfg.(names{ii})=overrides.(names{ii});
end
assert(cfg.water_depth_m==100 && cfg.source_clearance_m==1, ...
    'This entry fixes the newly approved physical installation.');
cfg.tx_depth_m=cfg.water_depth_m-cfg.source_clearance_m;
cfg.rx_range_m=(cfg.tx_depth_m-cfg.rx_depth_m)/tand(cfg.geometry_angle_deg);
assert(isfile(cfg.bellhop_exe),'Bellhop is missing.');
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
outputs=struct('received_png',fullfile(cfg.output_dir,'all_received_ray_paths.png'), ...
    'fan_png',fullfile(cfg.output_dir,'full_emission_ray_fan.png'), ...
    'data_mat',fullfile(cfg.output_dir,'source1m_plot_data.mat'));
if cfg.render_only
    saved=load(outputs.data_mat,'data','cfg');
    saved.cfg.show_figures=cfg.show_figures;
    local_plots(saved.data,saved.cfg,outputs);
    return;
end
profile_range=max(cfg.range_box_m,cfg.fan_range_m);
x=(-profile_range:cfg.surface_dx_m:profile_range).';
eta=0.045*cos(2*pi*x/6+0.30)+0.020*cos(2*pi*x/2.2-1.10) ...
    +0.010*cos(2*pi*x/0.9+2.00);
surface=table(x,eta,'VariableNames',{'x_m','surface_depth_m'});
writetable(surface,fullfile(cfg.output_dir,'fixed_surface_profile.csv'));
fprintf('A: 100 m bottom; Tx 99 m; Rx 3 m, range %.9g m.\n',cfg.rx_range_m);
base=fullfile(cfg.output_dir,'receiver_A');
local_run(cfg,x,eta,base,'A',[]);
a=local_arrivals([base '.arr']);
assert(abs(a.tx_depth_m-cfg.datum_shift_m-cfg.tx_depth_m)<1e-6 && ...
    abs(a.rx_depth_m-cfg.datum_shift_m-cfg.rx_depth_m)<1e-6 && ...
    abs(a.rx_range_m-cfg.rx_range_m)<1e-6,'Source/receiver location changed.');
arr=array2table(a.values,'VariableNames',{'amplitude','phase_deg','delay_s', ...
    'delay_imag_s','launch_angle_deg','receive_angle_deg','top_bounces','bottom_bounces'});
writetable(arr,fullfile(cfg.output_dir,'all_arrivals.csv'));
fprintf('E: collect every receiver-contributing record (not just strongest paths).\n');
base=fullfile(cfg.output_dir,'receiver_E');
local_run(cfg,x,eta,base,'E',[]);
[eigen,~]=local_rays([base '.ray']);
assert(~isempty(eigen),'No receiver-contributing trajectories.');
% Rounding printed E angles must not change the traced launch grid.
launch_grid=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),cfg.beam_count);
grid_index=round(([eigen.angle_deg]-launch_grid(1))/(launch_grid(2)-launch_grid(1)))+1;
assert(all(grid_index>=1 & grid_index<=cfg.beam_count),'E angle outside launch grid.');
assert(max(abs([eigen.angle_deg]-launch_grid(grid_index)))<2e-5,'E launch angle mismatch.');
indices=unique(grid_index); launch=launch_grid(indices);
fprintf('R: independently retrace %d launch angles for %d E records.\n',numel(launch),numel(eigen));
base=fullfile(cfg.output_dir,'receiver_retraced_R');
local_run(cfg,x,eta,base,'R',launch);
[rays,~]=local_rays([base '.ray']);
assert(numel(rays)==numel(launch),'R records missing.');
beam_rows=cell(numel(eigen),1); received=eigen;
for ii=1:numel(eigen)
    ix=find(indices==grid_index(ii),1);
    endpoint=eigen(ii).xy(:,end);
    [xy,miss]=local_prefix(rays(ix).xy,endpoint);
    % This is reconstruction quality, not an exact point-receiver hit test.
    assert(miss<1e-3,'E endpoint/R reconstruction residual %.6g m.',miss);
    received(ii).xy=xy-[0;cfg.datum_shift_m];
    rx=[cfg.rx_range_m;cfg.rx_depth_m];
    beam_rows{ii}=table(ii,eigen(ii).angle_deg,eigen(ii).top,eigen(ii).bottom, ...
        norm(received(ii).xy(:,end)-rx),local_distance(received(ii).xy,rx),miss, ...
        size(xy,2),'VariableNames',{'E_record','launch_angle_deg','top_bounces', ...
        'bottom_bounces','endpoint_to_rx_m','closest_to_rx_m','retrace_residual_m','vertices'});
end
beams=vertcat(beam_rows{:});
writetable(beams,fullfile(cfg.output_dir,'all_received_beams.csv'));
classes=unique([arr.top_bounces arr.bottom_bounces],'rows');
class_rows=cell(size(classes,1),1);
for ii=1:size(classes,1)
    top=classes(ii,1); bot=classes(ii,2);
    q=arr(arr.top_bounces==top & arr.bottom_bounces==bot,:);
    ne=sum(beams.top_bounces==top & beams.bottom_bounces==bot);
    assert(ne>0,'Arrival bounce class missing from E.');
    h=conj(sum(q.amplitude.*exp(1i*deg2rad(q.phase_deg) ...
        -1i*2*pi*cfg.frequency_hz*complex(q.delay_s,q.delay_imag_s))));
    class_rows{ii}=table(top,bot,height(q),ne,real(h),imag(h),abs(h), ...
        min(q.delay_s),max(q.delay_s),'VariableNames',{'top','bottom','arrivals', ...
        'E_records','H_real','H_imag','amplitude','min_delay_s','max_delay_s'});
end
class_table=vertcat(class_rows{:});
writetable(class_table,fullfile(cfg.output_dir,'path_classes.csv'));
fprintf('R fan: %d full trajectories, illustrative RBOX %.6g m.\n',cfg.fan_beams,cfg.fan_range_m);
base=fullfile(cfg.output_dir,'emission_fan_R');
fan_cfg=cfg; fan_cfg.range_box_m=cfg.fan_range_m;
local_run(fan_cfg,x,eta,base,'R',linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),cfg.fan_beams));
[fan,~]=local_rays([base '.ray']);
assert(numel(fan)==cfg.fan_beams,'Fan records missing.');
fan_rows=cell(numel(fan),1);
for ii=1:numel(fan)
    fan(ii).xy=fan(ii).xy-[0;cfg.datum_shift_m];
    xy=fan(ii).xy;
    assert(norm(xy(:,1)-[0;cfg.tx_depth_m])<1e-6,'Ray source moved.');
    assert(max(xy(2,:))<=cfg.water_depth_m+1e-4,'Ray crossed the physical bottom.');
    fan_rows{ii}=table(ii,fan(ii).angle_deg,fan(ii).top,fan(ii).bottom,size(xy,2), ...
        min(xy(1,:)),max(xy(1,:)),max(xy(2,:)),'VariableNames', ...
        {'ray','launch_angle_deg','top','bottom','vertices','min_x_m','max_x_m','max_z_m'});
end
fan_table=vertcat(fan_rows{:});
writetable(fan_table,fullfile(cfg.output_dir,'emission_fan_summary.csv'));
data=struct('surface',surface,'arrivals',arr,'received',received,'beams',beams, ...
    'classes',class_table,'fan',fan,'fan_table',fan_table);
save(outputs.data_mat,'data','cfg','outputs','-v7.3');
local_plots(data,cfg,outputs);
fprintf('DONE: %d arrivals, %d receiver E records, %d full fan rays.\n', ...
    height(arr),height(beams),height(fan_table));
disp(class_table);
end

function local_run(c,x,eta,base,mode,angles)
ns=c.beam_count; alphas=c.angle_limits_deg;
if ~isempty(angles), ns=numel(angles); alphas=angles; end
inside=abs(x)<=c.range_box_m+1e-8;
x=x(inside); eta=eta(inside);
bottom=c.water_depth_m+c.datum_shift_m; sspmax=bottom+0.1;
env=sprintf(['''100m water, physical source clearance 1m''\n%.17g\n1\n''CVW *''\n' ...
    '2 0 %.17g\n0 %.17g 0 1 0 0 /\n%.17g %.17g 0 1 0 0 /\n' ...
    '''A*'' 0\n%.17g %.17g 0 %.17g %.17g 0 /\n' ...
    '1\n%.17g /\n1\n%.17g /\n1\n%.17g /\n''%c  X''\n%d\n%s/\n%.17g %.17g %.17g\n'], ...
    c.frequency_hz,sspmax,c.c0_mps,sspmax,c.c0_mps,sspmax,c.bottom_c_mps, ...
    c.bottom_rho_g_cm3,c.bottom_alpha_db_lambda,c.tx_depth_m+c.datum_shift_m, ...
    c.rx_depth_m+c.datum_shift_m,c.rx_range_m/1000,mode,ns, ...
    sprintf('%.17g ',alphas),c.step_m,sspmax,c.range_box_m/1000);
ati=sprintf('''C''\n%d\n%s',numel(x),sprintf('%.17g %.17g\n',[x/1000 eta+c.datum_shift_m].'));
bty=sprintf('''L''\n2\n%.17g %.17g\n%.17g %.17g\n',-c.range_box_m/1000,bottom,c.range_box_m/1000,bottom);
ext='.ray'; if mode=='A', ext='.arr'; end
cached=isfile([base ext]) && isfile([base '.prt']) && contains(fileread([base '.prt']),'CPU Time');
for item={'.env',env;'.ati',ati;'.bty',bty}.'
    cached=cached && isfile([base item{1}]) && strcmp(fileread([base item{1}]),item{2});
end
if cached, return; end
local_write([base '.env'],env); local_write([base '.ati'],ati); local_write([base '.bty'],bty);
[folder,name]=fileparts(base); previous=pwd; cleanup=onCleanup(@()cd(previous));
cd(folder);
[status,message]=system(sprintf('"%s" "%s"',c.bellhop_exe,name));
local_write([base '.run.log'],message);
assert(status==0 && isfile([base ext]),'Bellhop failed: %s',message);
prt=fileread([base '.prt']);
assert(contains(prt,'CPU Time') && ~contains(prt,'FATAL ERROR'),'Incomplete Bellhop run.');
assert(~contains(prt,'Source below or too near') && ~contains(prt,'Source above or too near'), ...
    'Input preprocessor moved the physical source.');
end

function a=local_arrivals(file)
fid=fopen(file,'r'); assert(fid>=0,'Cannot open %s.',file);
cleanup=onCleanup(@()fclose(fid));
first=strtrim(fgetl(fid));
assert(~contains(first,'2D'),'This reader expects the saved legacy ASCII arrival format.');
v=sscanf(first,'%f');
assert(numel(v)==4 && all(v(2:4)==1),'Expected single-source/single-receiver arrival file.');
a=struct('frequency_hz',v(1),'tx_depth_m',fscanf(fid,'%f',1), ...
    'rx_depth_m',fscanf(fid,'%f',1),'rx_range_m',fscanf(fid,'%f',1),'values',[]);
fscanf(fid,'%d',1); n=fscanf(fid,'%d',1);
raw=fscanf(fid,'%f');
assert(numel(raw)==8*n,'Incomplete or unexpected arrival data.');
a.values=reshape(raw,8,n).';
end

function [rays,header]=local_rays(file)
fid=fopen(file,'r'); assert(fid>=0,'Cannot open %s.',file);
cleanup=onCleanup(@()fclose(fid));
fgetl(fid); frequency=fscanf(fid,'%f',1);
sources=fscanf(fid,'%d',3); beams=fscanf(fid,'%d',2);
top=fscanf(fid,'%f',1); bottom=fscanf(fid,'%f',1);
fgetl(fid); kind=strtrim(fgetl(fid));
assert(strcmp(kind,'''rz''') && all(sources==1),'Expected one-source 2-D rays.');
header=struct('frequency_hz',frequency,'beam_counts',beams, ...
    'top_depth_m',top,'bottom_depth_m',bottom);
rays=struct('angle_deg',{},'top',{},'bottom',{},'xy',{});
while true
    angle=fscanf(fid,'%f',1);
    if isempty(angle), break; end
    h=fscanf(fid,'%d',3);
    assert(numel(h)==3 && h(1)>=1,'Incomplete ray header.');
    xy=fscanf(fid,'%f',[2 h(1)]);
    assert(numel(xy)==2*h(1) && all(isfinite(xy),'all'),'Incomplete ray vertices.');
    rays(end+1)=struct('angle_deg',angle,'top',h(2),'bottom',h(3),'xy',xy); %#ok<AGROW>
end
% E header contains launched beam count, NOT number of written E records.
assert(isempty(strtrim(fscanf(fid,'%c'))),'Unexpected trailing ray data.');
end

function [prefix,miss]=local_prefix(xy,endpoint)
d=diff(xy,1,2); d2=sum(d.^2,1);
t=sum((endpoint-xy(:,1:end-1)).*d,1)./max(d2,realmin);
t=max(0,min(1,t)); q=xy(:,1:end-1)+d.*t;
[miss2,ix]=min(sum((q-endpoint).^2,1));
miss=sqrt(miss2); prefix=[xy(:,1:ix) q(:,ix)];
end

function distance=local_distance(xy,point)
[~,distance]=local_prefix(xy,point);
end

function [color,kind]=local_style(ray)
if ray.bottom>0
    kind=3; color=[0.53 0.26 0.68];
elseif ray.top>0
    kind=2; color=[0.9 0.36 0.1];
else
    kind=1; color=[0.05 0.40 0.72];
end
end


function local_plots(data,cfg,outputs)
local_one_plot(data.received,data,cfg,outputs.received_png,false);
local_one_plot(data.fan,data,cfg,outputs.fan_png,true);
end

function local_one_plot(rays,data,cfg,file,isfan)
visible='off'; if cfg.show_figures, visible='on'; end
fig=figure('Visible',visible,'Color','w','Position',[20 20 1850 970]);
cleanup=onCleanup(@()close(fig));
layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
minx=min(arrayfun(@(r)min(r.xy(1,:)),rays));
maxx=max(arrayfun(@(r)max(r.xy(1,:)),rays));
xl=[minx-0.02*max(1,maxx-minx) maxx+0.02*max(1,maxx-minx)];
% Draw all bottom classes first, then D/S so those overlaid lines stay visible.
order=[find([rays.bottom]>0) find([rays.bottom]==0 & [rays.top]>0) find([rays.top]==0 & [rays.bottom]==0)];
for pp=1:3
    ax=nexttile(layout); hold(ax,'on');
    for ii=order
        [color,~]=local_style(rays(ii));
        if isfan, color=0.6*color+0.4; width=0.35; else, width=0.55; end
        plot(ax,rays(ii).xy(1,:),rays(ii).xy(2,:),'Color',color,'LineWidth',width);
        if ~isfan
            plot(ax,rays(ii).xy(1,end),rays(ii).xy(2,end),'.','Color',color,'MarkerSize',3);
        end
    end
    plot(ax,data.surface.x_m,data.surface.surface_depth_m,'k-','LineWidth',1.3);
    yline(ax,cfg.water_depth_m,'Color',[0.35 0.24 0.10],'LineWidth',2);
    plot(ax,0,cfg.tx_depth_m,'kd','MarkerFaceColor',[1 .75 .1],'MarkerSize',8);
    plot(ax,cfg.rx_range_m,cfg.rx_depth_m,'kp','MarkerFaceColor',[.2 .85 .3],'MarkerSize',12);
    if pp==1
        xlim(ax,xl); ylim(ax,[-0.5 101.5]); title(ax,'Full stored paths (finite domain)');
        text(ax,xl(1)+.03*diff(xl),100.6,'Seabed z = 100 m','FontSize',10);
    elseif pp==2
        xlim(ax,[-2 max(3,cfg.rx_range_m+1)]); ylim(ax,[-.25 5]);
        title(ax,'Sea surface and receiver');
    else
        xlim(ax,[-2 2]); ylim(ax,[97.5 100.5]);
        title(ax,'Source: z = 99 m; clearance = 1 m');
        plot(ax,[-.5 -.5],[cfg.tx_depth_m cfg.water_depth_m],'k-','LineWidth',1.5);
        text(ax,-.45,99.5,'1 m','FontSize',12);
    end
    set(ax,'YDir','reverse','FontSize',11); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Horizontal x (m)'); ylabel(ax,'Physical depth z (m)');
end
if isfan
    head=sprintf('Complete launch fan: %d / %d rays; %.4f to %.1f deg', ...
        numel(rays),cfg.fan_beams,cfg.angle_limits_deg);
    sub=sprintf('4 kHz | water depth 100 m | illustrative RBOX +/- %.0f m | no first-surface clipping',cfg.fan_range_m);
else
    head=sprintf('All receiver-contributing paths: %d E records / %d A arrivals', ...
        numel(rays),height(data.arrivals));
    sub=sprintf('Geometry %.4g deg | Rx x=%.6g m | receiver RBOX +/- %.0f m | E beams are not exact point hits', ...
        cfg.geometry_angle_deg,cfg.rx_range_m,cfg.range_box_m);
end
title(layout,{head,sub, ...
    'Blue: no bounce | Orange: surface only | Purple: bottom involved | Diamond: Tx | Star: Rx | unequal axis scales'}, ...
    'FontSize',14,'FontWeight','bold');
exportgraphics(fig,file,'Resolution',150);
end

function local_write(file,content)
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0,'Cannot write %s.',file);
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end
