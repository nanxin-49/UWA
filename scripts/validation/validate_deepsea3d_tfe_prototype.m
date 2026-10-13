function state=validate_deepsea3d_tfe_prototype(stage,name)
%VALIDATE_DEEPSEA3D_TFE_PROTOTYPE Checkpointed independent TFE experiment.
arguments
    stage (1,1) string = "status"
    name (1,1) string = ""
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_tfe_prototype');file=fullfile(out,'prototype.mat');
if stage=="init"
    assert(~isfile(file),'TFE:Exists','Do not overwrite existing checkpoint.');
    a=load(fullfile(root,'results','validation','deepsea3d_pm_feasibility','feasibility.mat'),'state');
    b=load(fullfile(root,'results','validation','deepsea3d_g1_azimuth','azimuth.mat'),'state');
    assert(b.state.reference_ready && b.state.G1_status=="CONDITIONAL_G1_PASS");
    state=struct('pm_cfg',a.state.cfg,'frozen_sea',a.state.sea,'dual_cfg',b.state.cfg, ...
        'dual_reference',b.state.runs(b.state.reference_index).result,'runs',struct([]), ...
        'G1_status',"INCONCLUSIVE",'global_G1',"OPEN",'WP2_allowed',false,'FK_solves',0, ...
        'matlab_version',version);
else
    a=load(file,'state');state=a.state;
end
switch stage
    case "init"
        cfg=state.dual_cfg;cfg.profile.amplitudes_m(:)=0;
        opts=local_opts([16 .25 8 32 .8]);r=solve_deepsea3d_tfe(cfg,opts);
        bg=deepsea3d_flat(cfg,0);
        assert(abs(r.rough_correction)==0 && r.reflect_native==bg.reflect_native);
        state.flat=struct('result',r,'mirror_algebra_passed',true);
    case "run"
        [family,v,alpha]=local_job(name);
        assert(isempty(state.runs) || ~any(string({state.runs.name})==name), ...
            'TFE:Duplicate','Keep original result; choose a new named job.');
        if family=="dual",cfg=state.dual_cfg;else,cfg=state.pm_cfg;end
        cfg.profile.amplitudes_m=cfg.profile.amplitudes_m*alpha;
        [before,~]=memory;r=solve_deepsea3d_tfe(cfg,local_opts(v));[after,~]=memory;
        entry=struct('name',name,'family',family,'alpha',alpha,'settings',v,'result',r, ...
            'process_before_MiB',before.MemUsedMATLAB/2^20,'process_after_MiB',after.MemUsedMATLAB/2^20);
        if isempty(state.runs),state.runs=entry;else,state.runs(end+1)=entry;end
        fprintf('%s delta=(%.12g,%.12g) boundary %.3g half %.3g PDE %.3g halfPDE %.3g TBC %.3g wall %.3fs\n', ...
            name,real(r.rough_correction),imag(r.rough_correction),r.boundary_residual_relative, ...
            r.offgrid_boundary_residual_relative,r.PDE_residual_relative,r.offgrid_PDE_residual_relative, ...
            r.transparent_residual_relative,r.wall_s);
        if family=="dual"
            fprintf('Dual frozen-reference increment difference %.9g\n', ...
                abs(r.rough_correction-state.dual_reference.rough_correction)/abs(state.dual_reference.rough_correction));
        end
    case "status"
        disp(struct('runs',numel(state.runs),'G1',state.G1_status,'WP2',false));return
    otherwise
        error('TFE:Stage','Unknown stage.');
end
save(file,'state','-v7');local_csv(state,out);
end
function [family,v,alpha]=local_job(name)
alpha=1;family="dual";
switch name
    case "dual_base",v=[24 .25 12 32 .8];
    case "dual_order",v=[24 .25 16 32 .8];
    case "dual_order20",v=[24 .25 20 32 .8];
    case "dual_horizontal",v=[24 .125 16 32 .8];
    case "dual_vertical",v=[24 .25 16 40 .8];
    case "dual_top",v=[24 .25 16 40 1.2];
    case "dual_domain",v=[32 .25 16 32 .8];
    case "pm16",family="pm";v=[24 .25 16 32 .8];
    case "pm24",family="pm";v=[24 .25 24 32 .8];
    case "pm32",family="pm";v=[24 .25 32 32 .8];
    case "pm40",family="pm";v=[24 .25 40 32 .8];
    case "pm_scaled",family="pm";alpha=.25;v=[24 .25 16 32 .8];
    case "pm_horizontal",family="pm";v=[24 .125 32 32 .8];
    case "pm_vertical",family="pm";v=[24 .25 32 40 .8];
    case "pm_top",family="pm";v=[24 .25 32 40 1.2];
    case "pm_domain",family="pm";v=[32 .25 32 32 .8];
    otherwise,error('TFE:Job','Undefined bounded job.');
end
end
function opts=local_opts(v)
opts=struct('domain_halfwidth_m',v(1),'spacing_m',v(2),'order',v(3), ...
    'vertical_degree',v(4),'top_depth_m',v(5),'max_workspace_MiB',900);
end
function local_csv(s,out)
rows=struct([]);terms=struct([]);
for j=1:numel(s.runs)
    a=s.runs(j);r=a.result;
    row=struct('name',a.name,'family',a.family,'alpha',a.alpha,'L_m',r.domain_halfwidth_m, ...
        'h_m',r.spacing_m,'N',r.order,'M',r.vertical_degree,'b_m',r.top_depth_m, ...
        'delta_real',real(r.rough_correction),'delta_imag',imag(r.rough_correction), ...
        'boundary',r.boundary_residual_relative,'offgrid_boundary',r.offgrid_boundary_residual_relative, ...
        'PDE',r.PDE_residual_relative,'offgrid_PDE',r.offgrid_PDE_residual_relative, ...
        'TBC',r.transparent_residual_relative,'linear_residual',max(r.linear_residual), ...
        'cancellation',max(r.recursion_cancellation),'wall_s',r.wall_s,'cpu_s',r.cpu_s, ...
        'volume_MiB',r.volume_MiB,'workspace_estimate_MiB',r.estimated_workspace_MiB, ...
        'process_before_MiB',a.process_before_MiB,'process_after_MiB',a.process_after_MiB);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    for n=1:r.order
        row=struct('name',a.name,'n',n,'term_real',real(r.terms(n)),'term_imag',imag(r.terms(n)), ...
            'volume_norm',r.volume_coefficient_norms(n),'recursion_indicator',r.recursion_cancellation(n));
        if isempty(terms),terms=row;else,terms(end+1)=row;end %#ok<AGROW>
    end
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'runs.csv'));end
if ~isempty(terms),writetable(struct2table(terms),fullfile(out,'terms.csv'));end
end
