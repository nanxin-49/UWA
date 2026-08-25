function validation=validate_gaussian_sponge_vertical(overrides)
%VALIDATE_GAUSSIAN_SPONGE_VERTICAL Production Gaussian sponge engineering audit.

if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
gate_file=fullfile(cfg.window_gate_dir,'gaussian_window_convergence.mat');
if exist(gate_file,'file')~=2, error('Extended Gaussian window gate is missing.'); end
d=load(gate_file,'validation'); gate=d.validation;
if ~gate.hard_check_passed || ~gate.window_converged
    error('Gaussian PE-AS/window gate failed; sponge scan is blocked.');
end

Wcmp=cfg.comparison_reference_width_m; Ncmp=local_n(Wcmp,cfg.dx_target_m);
reference=vertical_channel_model(local_params(cfg,Wcmp,Ncmp,0,0.12,'slice'));
reference_field=reference.psifinal_xy;

count=numel(cfg.windows_m)*numel(cfg.sponge_ratios)*numel(cfg.alpha_values);
rows=repmat(local_row(),count,1); q=0; saved_fields=struct();
for ww=1:numel(cfg.windows_m)
    W=cfg.windows_m(ww); N=local_n(W,cfg.dx_target_m);
    base=vertical_channel_model(local_params(cfg,W,N,0,cfg.sponge_ratios(1),'slice'));
    for rr=1:numel(cfg.sponge_ratios)
        ratio=cfg.sponge_ratios(rr); e0=local_energy(base.psifinal_xy,base.x,base.y,W,ratio);
        for aa=1:numel(cfg.alpha_values)
            alpha=cfg.alpha_values(aa); q=q+1;
            if alpha==0, out=base; else, out=vertical_channel_model(local_params(cfg,W,N,alpha,ratio,'slice')); end
            e=local_energy(out.psifinal_xy,out.x,out.y,W,ratio);
            ref_on_grid=local_interp(reference_field,reference.x,reference.y,out.x,out.y);
            rows(q)=local_metrics(out.psifinal_xy,base.psifinal_xy,ref_on_grid, ...
                out.x,out.y,W,ratio,alpha,e,e0,cfg);
            if W==cfg.production_width_m && ratio==cfg.default_ratio && ...
                    ismember(alpha,[0 cfg.default_alpha])
                tag=sprintf('production_alpha_%g',alpha); tag=strrep(tag,'.','p');
                saved_fields.(tag)=struct('psi_xy',out.psifinal_xy,'x',out.x,'y',out.y);
            end
        end
    end
end
matrix=struct2table(rows);
targets=abs(matrix.axis_tl_change_db)<cfg.target_axis_tl_db & ...
    abs(matrix.axis_phase_change_rad)<cfg.target_axis_phase_rad & ...
    abs(matrix.center2_energy_change_db)<cfg.target_center_energy_db;
matrix.meets_center_targets=targets;
matrix.meets_edge_target=matrix.edge_energy_change_db<=-cfg.edge_suppression_target_db;
matrix.suitable_sponge=matrix.meets_center_targets & matrix.meets_edge_target & matrix.alpha_max_np_per_m>0;

production_default=matrix(matrix.width_m==cfg.production_width_m & ...
    matrix.sponge_ratio==cfg.default_ratio & matrix.alpha_max_np_per_m==cfg.default_alpha,:);
candidates=matrix(matrix.suitable_sponge,:);
recommendation=struct('use_sponge',false,'width_m',gate.reference_width_m, ...
    'sponge_ratio',NaN,'alpha_max_np_per_m',0,'reason','No nonzero sponge met all preregistered center and edge targets; use converged no-sponge window.');
if ~isempty(candidates)
    candidates.rank_score=abs(candidates.axis_tl_change_db)/cfg.target_axis_tl_db + ...
        abs(candidates.axis_phase_change_rad)/cfg.target_axis_phase_rad + ...
        abs(candidates.center2_energy_change_db)/cfg.target_center_energy_db + ...
        candidates.alpha_max_np_per_m;
    candidates=sortrows(candidates,{'width_m','rank_score','alpha_max_np_per_m'}, {'ascend','ascend','ascend'});
    best=candidates(1,:); recommendation=struct('use_sponge',true,'width_m',best.width_m, ...
        'sponge_ratio',best.sponge_ratio,'alpha_max_np_per_m',best.alpha_max_np_per_m, ...
        'reason','Weakest practical sponge meeting preregistered center and >=3 dB edge-suppression targets.');
end
validation=struct('schema_version','1.0.0','config',cfg,'window_gate',gate, ...
    'matrix',matrix,'production_default',production_default,'recommendation',recommendation, ...
    'reference_field',struct('width_m',Wcmp,'psi_xy',reference_field,'x',reference.x,'y',reference.y), ...
    'saved_fields',saved_fields);
validation.files=local_save(validation);
end

function cfg=local_defaults(root)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_gaussian_window_sponge'), ...
    'window_gate_dir',fullfile(root,'results','validation','pe_gaussian_window_sponge','extended_window'), ...
    'frequency_hz',4000,'c0_mps',1500,'sigma_src_m',0.3,'distance_m',97, ...
    'production_width_m',50,'windows_m',[50 80 128],'comparison_reference_width_m',160, ...
    'dx_target_m',50/256,'stepz_lamb',0.5, ...
    'sponge_ratios',[0.05 0.08 0.10 0.12 0.15 0.20], ...
    'alpha_values',[0 0.01 0.025 0.05 0.10 0.15 0.30], ...
    'default_ratio',0.12,'default_alpha',0.15,'target_axis_tl_db',0.1, ...
    'target_axis_phase_rad',0.02,'target_center_energy_db',0.1, ...
    'edge_suppression_target_db',3,'diagnostic_epsilon_db',1e-12);
end
function cfg=local_overrides(cfg,o), n=fieldnames(o); for i=1:numel(n), if ~isfield(cfg,n{i}), error('Unknown override: %s',n{i}); end, cfg.(n{i})=o.(n{i}); end, end
function N=local_n(W,dx), N=2*round(W/dx/2); end
function p=local_params(cfg,W,N,alpha,ratio,save_mode)
p=struct('f0',cfg.frequency_hz,'c0',cfg.c0_mps,'z_max',100,'z_tx',cfg.distance_m,'z_rx',0, ...
    'xw',W,'yw',W,'nx',N,'ny',N,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'stepz_lamb',cfg.stepz_lamb,'sigma_src_m',cfg.sigma_src_m,'source_mode','gaussian', ...
    'sponge_ratio',ratio,'alpha_max_np_per_m',alpha,'validation_allow_extended_sponge_ratio',true, ...
    'validation_allow_extended_window',W>100,'env_mode','uniform','enable_surface_reflection',false, ...
    'enable_bubbles',false,'doppler_fn',[],'enforce_1_over_R',false,'show_figures',false, ...
    'save_mode',save_mode,'nout',4,'use_gpu',false);
end
function e=local_energy(psi,x,y,W,ratio)
[X,Y]=meshgrid(x,y); p=abs(psi).^2; dA=abs(x(2)-x(1))*abs(y(2)-y(1)); rho=hypot(X,Y);
edge=abs(X)>=W*(0.5-ratio)|abs(Y)>=W*(0.5-ratio);
e=struct('center05',sum(p(rho<=0.5))*dA,'center1',sum(p(rho<=1))*dA, ...
    'center2',sum(p(rho<=2))*dA,'edge',sum(p(edge))*dA,'total',sum(p,'all')*dA);
end
function out=local_interp(psi,x,y,xq,yq), [Xq,Yq]=meshgrid(xq,yq); out=interp2(x,y,psi,Xq,Yq,'linear'); end
function r=local_metrics(sp,base,ref,x,y,W,ratio,alpha,e,e0,cfg)
[X,Y]=meshgrid(x,y); rho=hypot(X,Y); [~,ix]=min(abs(x)); [~,iy]=min(abs(y));
r=local_row(); r.width_m=W; r.nx=numel(x); r.dx_m=abs(x(2)-x(1)); r.sponge_ratio=ratio; r.alpha_max_np_per_m=alpha;
r.H_axis=sp(iy,ix); r.H0_axis=base(iy,ix); r.axis_amplitude_change_db=20*log10(abs(r.H_axis/r.H0_axis));
r.axis_tl_change_db=-r.axis_amplitude_change_db; r.axis_phase_change_rad=angle(r.H_axis*conj(r.H0_axis));
for radius=[0.5 1 2]
    m=rho<=radius; key=strrep(sprintf('center%.1f',radius),'.','p');
    r.([key '_complex_l2_vs_same_window'])=norm(sp(m)-base(m))/norm(base(m));
    r.([key '_complex_l2_vs_large_ref'])=norm(sp(m)-ref(m))/norm(ref(m));
end
r.center05_energy_change_db=10*log10(e.center05/e0.center05); r.center1_energy_change_db=10*log10(e.center1/e0.center1);
r.center2_energy_change_db=10*log10(e.center2/e0.center2); r.edge_energy_change_db=10*log10(e.edge/e0.edge);
r.total_energy_change_db=10*log10(e.total/e0.total); r.center2_energy=e.center2; r.edge_energy=e.edge; r.total_energy=e.total;
r.edge_center_gain=abs(r.edge_energy_change_db)/(abs(r.center2_energy_change_db)+cfg.diagnostic_epsilon_db);
r.axis_amplitude_error_vs_large_db=20*log10(abs(r.H_axis/ref(iy,ix))); r.axis_tl_error_vs_large_db=-r.axis_amplitude_error_vs_large_db;
r.axis_phase_error_vs_large_rad=angle(r.H_axis*conj(ref(iy,ix)));
end
function r=local_row()
r=struct('width_m',NaN,'nx',NaN,'dx_m',NaN,'sponge_ratio',NaN,'alpha_max_np_per_m',NaN, ...
    'H_axis',complex(NaN),'H0_axis',complex(NaN),'axis_amplitude_change_db',NaN,'axis_tl_change_db',NaN, ...
    'axis_phase_change_rad',NaN,'center0p5_complex_l2_vs_same_window',NaN,'center0p5_complex_l2_vs_large_ref',NaN, ...
    'center1p0_complex_l2_vs_same_window',NaN,'center1p0_complex_l2_vs_large_ref',NaN, ...
    'center2p0_complex_l2_vs_same_window',NaN,'center2p0_complex_l2_vs_large_ref',NaN, ...
    'center05_energy_change_db',NaN,'center1_energy_change_db',NaN,'center2_energy_change_db',NaN, ...
    'edge_energy_change_db',NaN,'total_energy_change_db',NaN,'center2_energy',NaN,'edge_energy',NaN,'total_energy',NaN, ...
    'edge_center_gain',NaN,'axis_amplitude_error_vs_large_db',NaN,'axis_tl_error_vs_large_db',NaN,'axis_phase_error_vs_large_rad',NaN);
end
function files=local_save(v)
out=v.config.output_dir; csv=fullfile(out,'gaussian_sponge_matrix.csv'); mat=fullfile(out,'gaussian_sponge_validation.mat');
writetable(v.matrix,csv); validation=v; save(mat,'validation','-v7.3'); files=struct('csv',csv,'mat',mat);
end
