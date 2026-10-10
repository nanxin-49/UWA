function output = validate_bellhop_fig12b(stage, overrides)
%VALIDATE_BELLHOP_FIG12B Native AT2020 deep-ocean signed-range audit.
% setup_vertical_project; validate_bellhop_fig12b('init');
% Repeated 'next' runs one checkpoint; 'render' plots saved numeric data.
% 'case' runs one explicit diagnostic override, without changing the plan.
% 'full' runs/reuses all 34 final cases, then regenerates figures and tables.
arguments
    stage (1,:) char = 'status'
    overrides (1,1) struct = struct()
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
cfg = struct('toolbox_root','E:/MISC/BELLHOP/AcousticsToolbox_2020', ...
    'output_dir',fullfile(root,'results','validation','bellhop_fig12b'), ...
    'tag','diagnostic','beam_count',50,'angles_deg',[-89.999 -89.94], ...
    'run_type','CG R','step_m',10,'range_m',linspace(-5,5,401), ...
    'depth_m',linspace(0,5000,1001),'water_depth_m',5000, ...
    'source_depth_m',4950,'frequency_hz',10000,'rbox_m',6, ...
    'munk_c0_mps',1500,'munk_axis_m',1300,'munk_scale_m',1300, ...
    'munk_epsilon',0.00737,'ssp_spacing_m',10, ...
    'bottom_cp_mps',1600,'bottom_rho_gcm3',1.8, ...
    'bottom_atten_db_wavelength',0.8,'uniform_control',false, ...
    'matched_control',false,'cerveny_type','MS', ...
    'epsilon_multiplier',1,'loop_range_m',5000,'image_count',1,'beam_window',5);
names=fieldnames(overrides);
for ii=1:numel(names)
    assert(isfield(cfg,names{ii}),'Unknown override %s',names{ii});
    cfg.(names{ii})=overrides.(names{ii});
end
assert(all(diff(cfg.range_m)>0) && numel(cfg.range_m)>1);
assert(max(abs(diff(cfg.range_m)-mean(diff(cfg.range_m))))<1e-9, ...
    'Native influence requires a regular range grid.');
assert(cfg.step_m>0 && cfg.beam_count>=1 && numel(cfg.angles_deg)==2);
assert(all(cfg.depth_m>=0 & cfg.depth_m<=cfg.water_depth_m));
addpath(fullfile(cfg.toolbox_root,'Matlab','ReadWrite'));
cfg.exe=fullfile(cfg.toolbox_root,'windows-bin-20201102','bellhop.exe');
assert(isfile(cfg.exe),'Native Bellhop executable missing.');
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
plan_file=fullfile(cfg.output_dir,'plan.mat');
switch lower(stage)
    case {'full','extend'}
        plan=local_plan(cfg); save(plan_file,'plan','cfg');
        if strcmpi(stage,'extend'), output=plan; return; end
        for ii=1:numel(plan)
            base=cfg; fields=fieldnames(plan(ii));
            for jj=1:numel(fields), base.(fields{jj})=plan(ii).(fields{jj}); end
            checkpoint=fullfile(cfg.output_dir,[plan(ii).tag '.mat']);
            if ~isfile(checkpoint)
                local_run(base);
            else
                local_check_checkpoint(checkpoint,base);
            end
        end
        output=generate_bellhop_fig12b_visuals(cfg.output_dir);
        output.audit=audit_bellhop_fig12b_artifacts(cfg.output_dir);
    case 'init'
        assert(~isfile(plan_file),'Existing plan retained; use next/status.');
        plan=local_plan(cfg);
        save(plan_file,'plan','cfg');
        output=plan; fprintf('Created %d reproducible cases.\n',numel(plan));
    case 'next'
        saved=load(plan_file,'plan','cfg');
        plan=saved.plan; base=saved.cfg;
        for ii=1:numel(plan)
            target=fullfile(base.output_dir,[plan(ii).tag '.mat']);
            if ~isfile(target)
                fields=fieldnames(plan(ii));
                for jj=1:numel(fields), base.(fields{jj})=plan(ii).(fields{jj}); end
                output=local_run(base); return
            end
        end
        output='All planned runs complete'; disp(output);
    case 'case'
        output=local_run(cfg);
    case 'render'
        output=generate_bellhop_fig12b_visuals(cfg.output_dir);
    case 'status'
        saved=load(plan_file,'plan'); plan=saved.plan;
        for ii=1:numel(plan)
            fprintf('%-30s %d\n',plan(ii).tag,isfile(fullfile(cfg.output_dir,[plan(ii).tag '.mat'])));
        end
        output=plan;
    otherwise
        error('Unknown stage %s',stage);
end
end

function plan=local_plan(cfg)
template=struct('tag','','beam_count',50,'angles_deg',cfg.angles_deg, ...
    'run_type','CG R','step_m',10,'depth_m',cfg.depth_m, ...
    'matched_control',false);
plan=repmat(template,1,12); k=0;
fans={[-89.999 -89.94],[-90.06 -89.94],[-179.999 179.999]};
labels={'initial','bilateral','wide'};
for ff=1:3
    for n=[50 500 2000]
        k=k+1; plan(k)=template;
        plan(k).tag=sprintf('%s_n%d_h10',labels{ff},n);
        plan(k).beam_count=n; plan(k).angles_deg=fans{ff};
    end
end
for h=[5 2.5 1.25]
    k=k+1; plan(k)=template;
    plan(k).tag=sprintf('wide_n2000_h%g',h);
    plan(k).beam_count=2000; plan(k).angles_deg=fans{3}; plan(k).step_m=h;
end
for h=[.625 .3125]
    plan(end+1)=local_case(template,sprintf('wide_n2000_h%g',h), ...
        struct('beam_count',2000,'angles_deg',fans{3},'step_m',h)); %#ok<AGROW>
end
plan(end+1)=local_case(template,'wide_n10000_h0.3125', ...
    struct('beam_count',10000,'angles_deg',fans{3},'step_m',.3125));
for n=[50 500 2000]
    delta=360/n;
    if mod(n,4)==0, angles=[-180+delta/2 180-delta/2];
    else, angles=[-180 180-delta]; end
    plan(end+1)=local_case(template,sprintf('sym_n%d_h0.3125',n), ...
        struct('beam_count',n,'angles_deg',angles,'step_m',.3125)); %#ok<AGROW>
end
plan(end+1)=local_case(template,'sym_n2000_h0.15625', ...
    struct('beam_count',2000,'angles_deg',[-179.91 179.91],'step_m',.15625));
for h=[1.25 .625 .3125]
    plan(end+1)=local_case(template,sprintf('bilateral_n2000_h%g',h), ...
        struct('beam_count',2000,'angles_deg',fans{2},'step_m',h)); %#ok<AGROW>
end
for n=[500 8000]
    plan(end+1)=local_case(template,sprintf('bilateral_n%d_h0.625',n), ...
        struct('beam_count',n,'angles_deg',fans{2},'step_m',.625)); %#ok<AGROW>
end
plan(end+1)=local_case(template,'left_only',struct('angles_deg',[-90.06 -90.001]));
plan(end+1)=local_case(template,'vertical_included',struct('beam_count',51,'angles_deg',fans{2}));
plan(end+1)=local_case(template,'strict_vertical_field',struct('beam_count',1,'angles_deg',[-90 -90]));
plan(end+1)=local_case(template,'strict_vertical_ray',struct('beam_count',1,'angles_deg',[-90 -90],'run_type','R'));
plan(end+1)=local_case(template,'rays_initial',struct('run_type','R'));
plan(end+1)=local_case(template,'rays_bilateral',struct('run_type','R','angles_deg',fans{2}));
ctrl=struct('beam_count',2000,'angles_deg',[-179.91 179.91], ...
    'step_m',.3125,'matched_control',true,'depth_m',[1000 3000 4900]);
plan(end+1)=local_case(template,'control_cart_2000',ctrl);
ctrl.run_type='Cg R'; plan(end+1)=local_case(template,'control_raycen_2000',ctrl);
ctrl.run_type='CG R'; ctrl.angles_deg=fans{2};
plan(end+1)=local_case(template,'control_narrow_2000',ctrl);
plan(end+1)=local_case(template,'surface_patch',struct('beam_count',2000, ...
    'angles_deg',fans{2},'step_m',.625,'depth_m',0:.01:1));
end

function p=local_case(template,tag,overrides)
p=template; p.tag=tag; fields=fieldnames(overrides);
for ii=1:numel(fields), p.(fields{ii})=overrides.(fields{ii}); end
end

function local_check_checkpoint(file,cfg)
s=load(file,'output'); old=s.output.cfg; fields=fieldnames(cfg);
for ii=1:numel(fields)
    name=fields{ii};
    % Early checkpoints predate this additional control flag (default false).
    if strcmp(name,'matched_control') && ~isfield(old,name), old.(name)=false; end
    % Gaussian options do not alter earlier geometric-beam environments.
    if any(strcmp(name,{'cerveny_type','epsilon_multiplier','loop_range_m','image_count','beam_window'})) ...
            && ~isfield(old,name) && (numel(old.run_type)<2 || ~any(old.run_type(2)=='CR'))
        old.(name)=cfg.(name);
    end
    assert(isfield(old,name) && isequaln(old.(name),cfg.(name)), ...
        'Checkpoint parameter %s differs: use a new output_dir (%s).',name,file);
end
end

function output=local_run(cfg)
target=fullfile(cfg.output_dir,[cfg.tag '.mat']);
assert(~isfile(target),'Checkpoint exists: %s',target);
case_dir=fullfile(cfg.output_dir,'bellhop');
if ~isfolder(case_dir), mkdir(case_dir); end
case_root=fullfile(case_dir,cfg.tag);
env_file=local_env(case_root,cfg);
% AT2020 FileRoot is CHARACTER(LEN=80). Set only the child shell directory;
% MATLAB's shared current directory is unchanged and all outputs are explicit.
assert(~contains(cfg.tag,{'"','&','|','<','>','%'}));
command=sprintf('cd /d "%s" && "%s" "%s"',case_dir,cfg.exe,cfg.tag);
t=tic; [status,command_output]=system(command); elapsed_s=toc(t);
fid=fopen([case_root '.console.txt'],'w');
assert(fid>=0); fprintf(fid,'%s',command_output); fclose(fid);
prt_file=[case_root '.prt'];
assert(isfile(prt_file),'Bellhop did not produce a print file.');
prt=fileread(prt_file);
assert(status==0 && ~contains(prt,'FATAL') && ~contains(prt,'Fatal'), ...
    'Native run failed: %s',command_output);
output=struct('cfg',cfg,'elapsed_s',elapsed_s,'exit_status',status, ...
    'command',command,'command_output',command_output, ...
    'env_file',env_file,'prt_file',prt_file, ...
    'storage_warnings',count(prt,'Insufficient storage'), ...
    'unbounded_warnings',count(prt,'Unbounded beam'), ...
    'matlab_version',version,'pressure',[],'range_m',[],'depth_m',[]);
if cfg.run_type(1)=='R'
    output.rays=local_rays([case_root '.ray']);
else
    [~,~,~,~,~,pos,p]=read_shd([case_root '.shd']);
    output.pressure=double(reshape(p,numel(pos.r.z),numel(pos.r.r)));
    output.range_m=double(pos.r.r(:).');
    output.depth_m=double(pos.r.z(:));
    assert(max(abs(output.range_m-cfg.range_m))<1e-5);
    assert(max(abs(output.depth_m-cfg.depth_m(:)))<1e-3);
    output.nonfinite_count=nnz(~isfinite(output.pressure));
    output.zero_fraction=nnz(output.pressure==0)/numel(output.pressure);
    zero_col=output.range_m==0;
    output.zero_axis_all_zero=any(zero_col) && all(output.pressure(:,zero_col)==0,'all');
    fprintf('%s: %.2fs, zero=%.5f, nonfinite=%d, storage=%d\n', ...
        cfg.tag,elapsed_s,output.zero_fraction,output.nonfinite_count,output.storage_warnings);
end
save(target,'output','-v7');
end

function env_file=local_env(case_root,cfg)
env_file=[case_root '.env'];
fid=fopen(env_file,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
z=(0:cfg.ssp_spacing_m:cfg.water_depth_m).';
eta=2*(z-cfg.munk_axis_m)/cfg.munk_scale_m;
c=cfg.munk_c0_mps*(1+cfg.munk_epsilon*(eta-1+exp(-eta)));
if cfg.uniform_control, c(:)=cfg.munk_c0_mps; end
fprintf(fid,'''Fig12b native audit - %s - assumed Munk/halfspace''\n',cfg.tag);
if cfg.matched_control
    z=[-10000;10000]; c=cfg.munk_c0_mps*ones(size(z));
    fprintf(fid,'%.12g\n1\n''CAW''\n',cfg.frequency_hz);
    fprintf(fid,'-10000 %.12g 0 1 0 /\n',cfg.munk_c0_mps);
    fprintf(fid,'2 0 10000\n');
else
    fprintf(fid,'%.12g\n1\n''SVW''\n%d 0 %.12g\n',cfg.frequency_hz,numel(z),cfg.water_depth_m);
end
fprintf(fid,'%.12g %.12g /\n',[z c].');
if cfg.matched_control
    fprintf(fid,'''A'' 0\n10000 %.12g 0 1 0 /\n',cfg.munk_c0_mps);
else
    fprintf(fid,'''A'' 0\n%.12g %.12g 0 %.12g %.12g /\n', ...
        cfg.water_depth_m,cfg.bottom_cp_mps,cfg.bottom_rho_gcm3,cfg.bottom_atten_db_wavelength);
end
fprintf(fid,'1\n%.12g /\n%d\n',cfg.source_depth_m,numel(cfg.depth_m));
fprintf(fid,'%.12g ',cfg.depth_m); fprintf(fid,'/\n%d\n',numel(cfg.range_m));
fprintf(fid,'%.12g ',cfg.range_m/1000); fprintf(fid,'/\n');
fprintf(fid,'''%s''\n%d\n%.12g %.12g /\n%.12g %.12g %.12g\n', ...
    cfg.run_type,cfg.beam_count,cfg.angles_deg,cfg.step_m,cfg.water_depth_m+100,cfg.rbox_m/1000);
if numel(cfg.run_type)>=2 && any(cfg.run_type(2)=='CR')
    % Native extra records; one true image avoids an assumed rigid-bottom image.
    assert(any(strcmp(cfg.cerveny_type,{'MS','MZ','FS','FZ'})));
    assert(cfg.epsilon_multiplier>0 && cfg.loop_range_m>0);
    fprintf(fid,'''%s'' %.12g %.12g\n%d %d ''P''\n',cfg.cerveny_type, ...
        cfg.epsilon_multiplier,cfg.loop_range_m/1000,cfg.image_count,cfg.beam_window);
end
clear cleanup
end

function rays=local_rays(file)
fid=fopen(file,'r'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fgetl(fid); fscanf(fid,'%f',1); dims=fscanf(fid,'%d',3);
assert(all(dims==1),'Ray reader expects one source.');
fscanf(fid,'%d',2); % number of elevation and bearing beams
fscanf(fid,'%f',2); fgetl(fid); fgetl(fid);
rays=struct('angle_deg',{},'top_bounces',{},'bottom_bounces',{},'rz_m',{});
while ~feof(fid)
    a=fscanf(fid,'%f',1); if isempty(a), break; end
    n=fscanf(fid,'%d',3); assert(numel(n)==3);
    rz=fscanf(fid,'%f',[2 n(1)]).'; assert(size(rz,1)==n(1));
    rays(end+1)=struct('angle_deg',a,'top_bounces',n(2), ...
        'bottom_bounces',n(3),'rz_m',rz); %#ok<AGROW>
end
clear cleanup
end
