function receipt=render_r2_three_seed_surfaces(output_dir)
%RENDER_R2_THREE_SEED_SURFACES Frozen Fourier/C2 geometry only, no acoustics.
arguments
    output_dir (1,:) char
end
assert(~isempty(regexp(output_dir,'^([A-Za-z]:[\\/]|/)','once')),'Use an explicit absolute output directory');
assert(~isfolder(output_dir)); mkdir(output_dir);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldpath=path; pg=onCleanup(@()path(oldpath)); addpath(fullfile(root,'scripts','validation','support'));
timer=tic;
case_root=fullfile(root,'results','comparisons','piston_2d_pm','r2_saved_closeout_20261117_01');
master=fullfile(case_root,'committed.json');
assert(strcmp(r1_bie_continuation_io('sha',master),'439d497821065817a2cd1cd009ee20ad239e5dcb60c7b861864b3c873fc63a24'));
seal=jsondecode(fileread(master)); pre=jsondecode(fileread(fullfile(case_root,'source_preflight.json')));
local_bound(case_root,seal.files,'source_preflight.json');
x=(-6000:6000)'*.01; depth=zeros(numel(x),3); elevation=depth;
profile=cell(1,3); sources=struct('seed',{},'path',{},'sha256',{},'saved_geometry_difference_m',{}); stats=table;
for j=1:3
    seed=260000+j;
    if j==1
        rec=pre.sources([pre.sources.seed]==seed); file=fullfile(rec.work,'common.mat');
        list=dir(fullfile(rec.source,'retrieval','SHA256SUMS_*.txt')); assert(isscalar(list));
        listfile=fullfile(list.folder,list.name); assert(strcmp(r1_bie_continuation_io('sha',listfile),rec.retrieval_manifest_sha256));
        rows=splitlines(strtrim(fileread(listfile))); found=false;
        for k=1:numel(rows)
            tok=regexp(char(rows(k)),'^([0-9a-f]{64})  (.+)$','tokens','once'); assert(~isempty(tok));
            rel=regexprep(tok{2},'^\./','');
            if strcmp(rel,'finalization/attempts/attempt_0001/work/common.mat')
                assert(strcmp(r1_bie_continuation_io('sha',file),tok{1})); found=true;
            end
        end
        assert(found);
    else
        relative=sprintf('seed%d/work/common.mat',seed); local_bound(case_root,seal.files,relative);
        file=fullfile(case_root,strrep(relative,'/',filesep));
    end
    q=load(file,'cfg','profile','sx','eta'); profile{j}=q.profile;
    assert(q.cfg.seed==seed && height(q.profile)==12 && q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000);
    [saved,~,~]=closeout_s3a_pm_surface(q.profile,q.sx);
    difference=max(abs(saved-q.eta)); assert(difference<1e-12);
    depth(:,j)=closeout_s3a_pm_surface(q.profile,x); elevation(:,j)=-depth(:,j);
    assert(all(depth(abs(x)>=50,j)==0));
    sources(j)=struct('seed',seed,'path',file,'sha256',r1_bie_continuation_io('sha',file),'saved_geometry_difference_m',difference);
    stats=[stats;table(seed,min(elevation(:,j)),max(elevation(:,j)),sqrt(mean(elevation(abs(x)<=40.0390625,j).^2)), ...
        'VariableNames',{'seed','sampled_min_elevation_m','sampled_max_elevation_m','RMS_over_receiver_span_m'})]; %#ok<AGROW>
end
assert(isequal(profile{1}.k_rad_per_m,profile{2}.k_rad_per_m,profile{3}.k_rad_per_m));
assert(isequal(profile{1}.Sk_m3,profile{2}.Sk_m3,profile{3}.Sk_m3));
T=table(x,depth(:,1),depth(:,2),depth(:,3),elevation(:,1),elevation(:,2),elevation(:,3), ...
    'VariableNames',{'x_m','depth_260001_m','depth_260002_m','depth_260003_m', ...
    'elevation_260001_m','elevation_260002_m','elevation_260003_m'});
writetable(T,fullfile(output_dir,'surface_profiles.csv')); writetable(stats,fullfile(output_dir,'surface_geometry_stats.csv'));
save(fullfile(output_dir,'surface_profiles.mat'),'x','depth','elevation','profile','sources','stats','-v7');
f=figure('Visible','off','Units','pixels','Position',[40 40 1800 1400],'Resize','off','Color','w','MenuBar','none','ToolBar','none');
fg=onCleanup(@()local_close(f)); colors=[0 .447 .698;.009 .55 .35;.835 .369 0]; styles={'-','--','-.'};
annotation(f,'textbox',[.065 .960 .9 .028],'String','R1 + R2 | Three frozen PM sea surfaces', ...
    'EdgeColor','none','FontName','Arial','FontSize',16,'HorizontalAlignment','center','Interpreter','none');
annotation(f,'textbox',[.065 .935 .9 .025],'String','Elevation eta = -z_s; original coefficients, 160 m period, 12 modes; vertical scale enlarged for visibility', ...
    'EdgeColor','none','FontName','Arial','FontSize',11,'HorizontalAlignment','center','Interpreter','none');
positions=[150 1030 1330 255;150 710 1330 255;150 390 1330 255;150 70 1330 255];
limit=ceil(max(abs(elevation),[],'all')*10)/10+.05;
for panel=1:4
    ax=axes('Parent',f,'Units','pixels','Position',positions(panel,:),'PositionConstraint','innerposition', ...
        'LooseInset',[0 0 0 0],'FontName','Arial','FontSize',11,'Box','on'); hold(ax,'on');
    for range=[-50 -42;42 50].'
        patch(ax,[range(1) range(2) range(2) range(1)],[-limit -limit limit limit], ...
            [.93 .94 .94],'EdgeColor','none','HandleVisibility','off');
    end
    yline(ax,0,':','Color',[.5 .5 .5],'HandleVisibility','off');
    which_seeds=1:3; if panel>1, which_seeds=panel-1; end
    for j=which_seeds
        line=plot(ax,x,elevation(:,j),styles{j},'Color',colors(j,:),'LineWidth',1.3,'DisplayName',sprintf('seed %d',260000+j));
        line.UserData=struct('role','surface_profile','seed',260000+j,'panel',panel);
    end
    ax.XLim=[-60 60];ax.YLim=[-limit limit];ax.XTick=-60:10:60;ax.YTick=linspace(-limit,limit,5);
    ax.XTickLabel=arrayfun(@(v)sprintf('%g',v),ax.XTick,'UniformOutput',false);
    ax.YTickLabel=arrayfun(@(v)sprintf('%.3g',v),ax.YTick,'UniformOutput',false);
    ax.XAxis.Exponent=0;ax.YAxis.Exponent=0;grid(ax,'on');
    ylabel(ax,'Surface elevation eta (m)','Interpreter','none');xlabel(ax,'Horizontal coordinate x (m)');
    lg=legend(ax,'Location','northeast','AutoUpdate','off','FontSize',11);lg.Units='pixels';
    lg.Position=[1500 positions(panel,2)+163 280 82];
end
annotation(f,'textbox',[.835 .745 .15 .105],'String',{'Grey: C2 transition', '|x| <= 42 m: unchanged', '|x| >= 50 m: flat', '', 'Receiver span:', '-40.039 to 40.039 m', '', 'Geometry only; no solve'}, ...
    'EdgeColor','none','FontName','Arial','FontSize',10,'Interpreter','none');
drawnow; stem='F06_three_seed_surfaces'; savefig(f,fullfile(output_dir,[stem '.fig']));
exportgraphics(f,fullfile(output_dir,[stem '.png']),'Resolution',160,'Width',1800,'Height',1400,'Units','pixels'); clear fg;
% Independent Fourier sum and C2 reconstruction, including original saved samples.
independent=zeros(size(elevation)); r=abs(x); window=double(r<=42); ii=r>42 & r<50;t=(r(ii)-42)/8;
window(ii)=1-10*t.^3+15*t.^4-6*t.^5;
for j=1:3
    p=profile{j}; raw=zeros(size(x));
    for k=1:height(p), raw=raw+p.cos_coeff_m(k)*cos(p.k_rad_per_m(k)*x)+p.sin_coeff_m(k)*sin(p.k_rad_per_m(k)*x); end
    independent(:,j)=-raw.*window;
end
assert(max(abs(independent-elevation),[],'all')<1e-12);
f=openfig(fullfile(output_dir,[stem '.fig']),'new','invisible');fg=onCleanup(@()local_close(f)); count=0;
lines=findall(f,'Type','line');
for j=1:numel(lines)
    u=lines(j).UserData;if ~isstruct(u)||~isfield(u,'role')||~strcmp(u.role,'surface_profile'),continue;end
    index=u.seed-260000; assert(isequal(lines(j).XData(:),x)&&isequal(lines(j).YData(:),elevation(:,index)));count=count+1;
end
assert(count==6);original=imread(fullfile(output_dir,[stem '.png']));pixels=zeros(1,2);
for pass=1:2
    if pass==2,clear fg;f=openfig(fullfile(output_dir,[stem '.fig']),'new','invisible');fg=onCleanup(@()local_close(f));end
    f.Units='pixels';f.Position=[40 40 1800 1400];drawnow;
    name=sprintf('%s_reopen%d',stem,pass);savefig(f,fullfile(output_dir,[name '.fig']));
    exportgraphics(f,fullfile(output_dir,[name '.png']),'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
    audit=imread(fullfile(output_dir,[name '.png']));assert(isequal(size(original),size(audit)));
    pixels(pass)=nnz(original~=audit);assert(pixels(pass)==0);
end
receipt=struct('status','SAVED_GEOMETRY_AND_FIG_PNG_PASS','failed_checks',0,'sources',sources,'MATLAB_version',version, ...
    'sample_spacing_m',.01,'display_interval_m',[-60 60],'checked_FIG_profiles',count,'PNG_channel_differences',pixels, ...
    'independent_fourier_difference_m',max(abs(independent-elevation),[],'all'),'elapsed_s',toc(timer), ...
    'new_acoustic_solves',0,'SSH_connections',0,'representation','elevation eta=-physical depth z_s; no amplitude normalization or new random samples');
r1_bie_manual_io('write',fullfile(output_dir,'receipt.json'),receipt);
fprintf('THREE_SEED_SURFACES status=%s elapsed_s=%.6f\n',receipt.status,receipt.elapsed_s);
end
function local_bound(root,files,rel)
ix=find(strcmp({files.path},rel));assert(isscalar(ix));
assert(strcmp(r1_bie_continuation_io('sha',fullfile(root,strrep(rel,'/',filesep))),files(ix).sha256));
end
function local_close(f)
if isgraphics(f),close(f);end
end
