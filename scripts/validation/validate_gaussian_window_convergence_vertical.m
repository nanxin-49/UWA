function validation=validate_gaussian_window_convergence_vertical(overrides)
%VALIDATE_GAUSSIAN_WINDOW_CONVERGENCE_VERTICAL Production Gaussian PE/AS/window audit.

if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

hard_rows=repmat(local_hard_row(),numel(cfg.distances_m),1);
for ii=1:numel(cfg.distances_m)
    L=cfg.distances_m(ii); W=cfg.hard_check_width_m; N=local_n(W,cfg.dx_target_m);
    out=vertical_channel_model(local_params(cfg,W,N,L,0,0.12,'slice',4));
    [X,Y]=meshgrid(out.x,out.y);
    [psi0,source_meta]=gaussian_source_initial_field_vertical(X,Y,out.config);
    psi_as=exact_angular_spectrum_one_step_vertical(psi0,W,W,cfg.frequency_hz,cfg.c0_mps,L);
    hard_rows(ii)=local_hard_metrics(out.psifinal_xy,psi_as,out.x,out.y,L,cfg.center_radii_m);
    hard_rows(ii).source_meta={source_meta};
end
hard_table=struct2table(hard_rows);
hard_pass=all(hard_table.full_complex_l2<=cfg.pe_as_full_l2_limit & ...
    hard_table.center_complex_l2_max<=cfg.pe_as_center_l2_limit & ...
    abs(hard_table.axis_amplitude_error_db)<=cfg.pe_as_axis_amplitude_limit_db & ...
    abs(hard_table.axis_phase_error_rad)<=cfg.pe_as_axis_phase_limit_rad);
if ~hard_pass
    validation=struct('schema_version','1.0.0','config',cfg,'source_audit',source_meta, ...
        'hard_check_table',hard_table,'hard_check_passed',false,'window_converged',false);
    files=local_save(validation,table(),table(),root);
    validation.files=files; %#ok<STRNU>
    error('Gaussian PE-AS hard check failed; window and sponge work stopped.');
end

fields=cell(numel(cfg.distances_m),numel(cfg.widths_m));
for ll=1:numel(cfg.distances_m)
    L=cfg.distances_m(ll);
    for ww=1:numel(cfg.widths_m)
        W=cfg.widths_m(ww); N=local_n(W,cfg.dx_target_m);
        out=vertical_channel_model(local_params(cfg,W,N,L,0,0.12,'slice',4));
        fields{ll,ww}=struct('psi',out.psifinal_xy,'x',out.x,'y',out.y, ...
            'psiout',out.psiout,'z_track',out.z_track,'config',out.config);
    end
end

rows=repmat(local_window_row(),numel(cfg.distances_m)*numel(cfg.widths_m)*numel(cfg.center_radii_m),1); q=0;
profile_rows=repmat(local_profile_row(),numel(cfg.distances_m)*numel(cfg.widths_m)*5,1); qp=0;
for ll=1:numel(cfg.distances_m)
    ref=fields{ll,end};
    for ww=1:numel(cfg.widths_m)
        cur=fields{ll,ww};
        ref_on_cur=local_interp_field(ref,cur.x,cur.y);
        for rr=1:numel(cfg.center_radii_m)
            q=q+1; rows(q)=local_region_metrics(cur.psi,ref_on_cur,cur.x,cur.y, ...
                cfg.distances_m(ll),cfg.widths_m(ww),cfg.center_radii_m(rr));
        end
        planes=local_plane_metrics(cur,cfg.widths_m(ww),cfg.distances_m(ll));
        for pp=1:numel(planes)
            qp=qp+1; profile_rows(qp)=planes(pp);
        end
    end
end
window_table=struct2table(rows); propagation_table=struct2table(profile_rows(1:qp));

last_pair=window_table(window_table.width_m==cfg.widths_m(end-1),:);
window_converged=all(last_pair.complex_l2<=cfg.window_complex_l2_limit & ...
    last_pair.max_amplitude_error_db<=cfg.window_max_amplitude_limit_db & ...
    last_pair.phase_rms_rad<=cfg.window_phase_rms_limit_rad);
reference_width_m=NaN;
if window_converged
    for ww=1:numel(cfg.widths_m)
        tail=window_table(window_table.width_m>=cfg.widths_m(ww),:);
        if all(tail.complex_l2<=cfg.window_complex_l2_limit & ...
                tail.max_amplitude_error_db<=cfg.window_max_amplitude_limit_db & ...
                tail.phase_rms_rad<=cfg.window_phase_rms_limit_rad)
            reference_width_m=cfg.widths_m(ww); break
        end
    end
end
validation=struct('schema_version','1.0.0','config',cfg,'source_audit',source_meta, ...
    'hard_check_table',hard_table,'hard_check_passed',hard_pass, ...
    'window_table',window_table,'propagation_table',propagation_table, ...
    'window_converged',window_converged,'reference_width_m',reference_width_m, ...
    'reference_definition','largest tested no-sponge field; accepted only if penultimate window passes preregistered limits');
validation.files=local_save(validation,window_table,propagation_table,root);
if ~window_converged
    error('No-sponge Gaussian field did not converge at the largest tested windows; sponge optimization stopped.');
end
end

function cfg=local_defaults(root)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_gaussian_window_sponge'), ...
    'frequency_hz',4000,'c0_mps',1500,'sigma_src_m',0.3, ...
    'distances_m',[20 40 70 97 100],'widths_m',[16 24 32 48 50 64 80 100], ...
    'production_width_m',50,'hard_check_width_m',100,'dx_target_m',50/256, ...
    'center_radii_m',[0.5 1 2],'stepz_lamb',0.5, ...
    'pe_as_full_l2_limit',1e-10,'pe_as_center_l2_limit',1e-10, ...
    'pe_as_axis_amplitude_limit_db',1e-9,'pe_as_axis_phase_limit_rad',1e-10, ...
    'window_complex_l2_limit',1e-3,'window_max_amplitude_limit_db',0.05, ...
    'window_phase_rms_limit_rad',0.01);
end
function cfg=local_overrides(cfg,o), n=fieldnames(o); for i=1:numel(n), if ~isfield(cfg,n{i}), error('Unknown override: %s',n{i}); end, cfg.(n{i})=o.(n{i}); end, end
function N=local_n(W,dx), N=2*round(W/dx/2); end
function p=local_params(cfg,W,N,L,alpha,ratio,save_mode,nout)
p=struct('f0',cfg.frequency_hz,'c0',cfg.c0_mps,'z_max',max(100,L),'z_tx',L, ...
    'z_rx',0,'xw',W,'yw',W,'nx',N,'ny',N,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'stepz_lamb',cfg.stepz_lamb,'sigma_src_m',cfg.sigma_src_m,'source_mode','gaussian', ...
    'sponge_ratio',ratio,'alpha_max_np_per_m',alpha,'env_mode','uniform', ...
    'enable_surface_reflection',false,'enable_bubbles',false,'doppler_fn',[], ...
    'validation_allow_extended_window',W>100, ...
    'enforce_1_over_R',false,'show_figures',false,'save_mode',save_mode,'nout',nout,'use_gpu',false);
end
function r=local_hard_metrics(pe,as,x,y,L,radii)
[X,Y]=meshgrid(x,y); rho=hypot(X,Y); d=pe-as; [~,ix]=min(abs(x)); [~,iy]=min(abs(y));
c=zeros(size(radii)); for i=1:numel(radii), m=rho<=radii(i); c(i)=norm(d(m))/norm(as(m)); end
r=local_hard_row(); r.distance_m=L; r.full_complex_l2=norm(d(:))/norm(as(:));
r.center_complex_l2_max=max(c); r.axis_amplitude_error_db=20*log10(abs(pe(iy,ix)/as(iy,ix)));
r.axis_phase_error_rad=angle(pe(iy,ix)*conj(as(iy,ix)));
end
function r=local_region_metrics(cur,ref,x,y,L,W,radius)
[X,Y]=meshgrid(x,y); m=hypot(X,Y)<=radius; a=abs(cur(m)); b=abs(ref(m)); valid=b>max(b)*1e-12;
ratio=a(valid)./b(valid); phase=angle(cur(m).*conj(ref(m)));
r=local_window_row(); r.distance_m=L; r.width_m=W; r.nx=numel(x); r.dx_m=abs(x(2)-x(1)); r.radius_m=radius;
r.complex_l2=norm(cur(m)-ref(m))/norm(ref(m)); r.mean_amplitude_error_db=mean(abs(20*log10(ratio)));
r.max_amplitude_error_db=max(abs(20*log10(ratio))); r.phase_rms_rad=sqrt(mean(phase.^2));
r.center_energy=sum(abs(cur(m)).^2)*r.dx_m^2; r.reference_center_energy=sum(abs(ref(m)).^2)*r.dx_m^2;
r.center_energy_error_db=10*log10(r.center_energy/r.reference_center_energy);
[~,ix]=min(abs(x)); [~,iy]=min(abs(y)); r.axis_amplitude_error_db=20*log10(abs(cur(iy,ix)/ref(iy,ix)));
r.axis_tl_error_db=-r.axis_amplitude_error_db; r.axis_phase_error_rad=angle(cur(iy,ix)*conj(ref(iy,ix)));
end
function out=local_interp_field(ref,x,y)
[Xq,Yq]=meshgrid(x,y); out=interp2(ref.x,ref.y,ref.psi,Xq,Yq,'linear');
end
function rows=local_plane_metrics(c,W,L)
N=size(c.psiout,1); rows=repmat(local_profile_row(),N+1,1); dx=abs(c.x(2)-c.x(1));
[X,Y]=meshgrid(c.x,c.y); psi0=gaussian_source_initial_field_vertical(X,Y,c.config);
numstep=numel(c.z_track)-1; idx=[0 round(numstep*(1:N)/N)]; fields=cell(N+1,1); fields{1}=psi0;
for i=1:N, fields{i+1}=squeeze(c.psiout(i,:,:)); end
for i=1:N+1
    p=abs(fields{i}).^2; center=hypot(X,Y)<=2; edge5=abs(X)>=0.45*W|abs(Y)>=0.45*W; edge10=abs(X)>=0.40*W|abs(Y)>=0.40*W;
    boundary=abs(X)>=0.49*W|abs(Y)>=0.49*W; rows(i)=local_profile_row(); rows(i).distance_m=L;
    rows(i).width_m=W; rows(i).propagated_m=L*idx(i)/numstep; rows(i).center_energy=sum(p(center))*dx^2;
    rows(i).edge5_energy=sum(p(edge5))*dx^2; rows(i).edge10_energy=sum(p(edge10))*dx^2;
    rows(i).total_energy=sum(p,'all')*dx^2; rows(i).boundary_max_relative_db=20*log10(max(abs(fields{i}(boundary)))/max(abs(fields{i}(:))));
end
end
function r=local_hard_row(), r=struct('distance_m',NaN,'full_complex_l2',NaN,'center_complex_l2_max',NaN,'axis_amplitude_error_db',NaN,'axis_phase_error_rad',NaN,'source_meta',{{}}); end
function r=local_window_row(), r=struct('distance_m',NaN,'width_m',NaN,'nx',NaN,'dx_m',NaN,'radius_m',NaN,'complex_l2',NaN,'mean_amplitude_error_db',NaN,'max_amplitude_error_db',NaN,'phase_rms_rad',NaN,'center_energy',NaN,'reference_center_energy',NaN,'center_energy_error_db',NaN,'axis_amplitude_error_db',NaN,'axis_tl_error_db',NaN,'axis_phase_error_rad',NaN); end
function r=local_profile_row(), r=struct('distance_m',NaN,'width_m',NaN,'propagated_m',NaN,'center_energy',NaN,'edge5_energy',NaN,'edge10_energy',NaN,'total_energy',NaN,'boundary_max_relative_db',NaN); end
function files=local_save(v,w,p,root)
out=v.config.output_dir; mat=fullfile(out,'gaussian_window_convergence.mat'); hard=fullfile(out,'gaussian_pe_as_hard_check.csv'); wc=fullfile(out,'gaussian_window_convergence.csv'); pc=fullfile(out,'gaussian_propagation_energy.csv');
writetable(v.hard_check_table,hard); if ~isempty(w), writetable(w,wc); end; if ~isempty(p), writetable(p,pc); end
validation=v; save(mat,'validation','-v7.3'); files=struct('mat',mat,'hard_csv',hard,'window_csv',wc,'propagation_csv',pc,'report',fullfile(root,'reports','pe_gaussian_window_sponge_validation_report.md'));
end
