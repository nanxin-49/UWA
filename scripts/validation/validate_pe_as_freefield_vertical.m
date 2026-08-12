function validation = validate_pe_as_freefield_vertical(overrides)
%VALIDATE_PE_AS_FREEFIELD_VERTICAL Validate production PE in free space.
% Level 1 is a hard implementation check against an independent one-step
% angular spectrum. Level 2 compares source-referenced PE with a spherical
% wave without any fitted complex scale.

if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
run_full_matrix=logical(cfg.run_full_matrix);

grid_cases=cfg.grid_cases;
field_rows=repmat(local_field_row(),numel(grid_cases),1);
field_data=cell(numel(grid_cases),1);
for ii=1:numel(grid_cases)
    g=grid_cases(ii);
    p=local_params(cfg,cfg.f_ref_hz,cfg.representative_L_m, ...
        cfg.representative_s0_m,0,g.nx,g.width_m,g.stepz_lamb,'slice',0);
    out=vertical_channel_model(p);
    [X,Y]=meshgrid(out.x,out.y);
    [psi0,~]=virtual_point_source_initial_field_vertical(X,Y,cfg.f_ref_hz,out.config);
    [psi_as,as_meta]=exact_angular_spectrum_one_step_vertical(psi0,g.width_m, ...
        g.width_m,cfg.f_ref_hz,cfg.c0_mps,cfg.representative_L_m);
    pe=out.psifinal_xy; delta=pe-psi_as;
    rho=hypot(X,Y); trusted=rho<=cfg.trusted_radius_m;
    phase_mask=trusted & abs(psi_as)>=max(abs(psi_as(:)))*cfg.phase_floor_ratio;
    k=2*pi*cfg.f_ref_hz/cfg.c0_mps;
    R=hypot(cfg.representative_s0_m+cfg.representative_L_m,rho);
    h_source_pe=pe.*exp(1i*k*(cfg.representative_s0_m+cfg.representative_L_m));
    h_exact=cfg.green_amplitude.*exp(1i*k*R)./R;
    row=local_field_row(); row.case_name=string(g.name); row.nx=g.nx;
    row.width_m=g.width_m; row.stepz_lamb=g.stepz_lamb;
    row.pe_as_rel_l2=norm(delta(:))/norm(psi_as(:));
    row.pe_as_max_abs=max(abs(delta(:)));
    row.pe_as_trusted_rel_l2=norm(delta(trusted))/norm(psi_as(trusted));
    row.pe_as_phase_rms_rad=local_phase_rms(pe,psi_as,phase_mask);
    row.pe_exact_trusted_rel_l2=norm(h_source_pe(trusted)-h_exact(trusted))/norm(h_exact(trusted));
    row.pe_exact_phase_rms_rad=local_phase_rms(h_source_pe,h_exact,phase_mask);
    row.pe_exact_tl_rms_db=local_tl_rms(h_source_pe,h_exact,phase_mask);
    row.invariant_error=max(abs(out.H_f(:)-out.H_direct_f(:)-out.H_reflect_f(:)));
    row.max_reflection=max(abs(out.H_reflect_f(:)));
    field_rows(ii)=row;
    field_data{ii}=struct('x_m',out.x,'y_m',out.y,'psi_pe',pe,'psi_as',psi_as, ...
        'h_source_pe',h_source_pe,'h_exact',h_exact,'trusted_mask',trusted, ...
        'phase_mask',phase_mask,'as_meta',as_meta,'source_meta',out.source_meta);
end
field_table=struct2table(field_rows);

if run_full_matrix
axis_rows=repmat(local_axis_row(),numel(cfg.s0_values_m)*numel(cfg.L_values_m)*numel(cfg.frequencies_hz),1);
q=0;
for ss=1:numel(cfg.s0_values_m)
    for ll=1:numel(cfg.L_values_m)
        p=local_params(cfg,cfg.frequencies_hz,cfg.L_values_m(ll),cfg.s0_values_m(ss), ...
            0,cfg.axis_nx,cfg.axis_width_m,cfg.axis_stepz_lamb,'rx_only',cfg.axis_alpha_max_np_per_m);
        out=vertical_channel_model(p); f=out.f_axis(:); k=2*pi*f/cfg.c0_mps;
        h_plane=out.H_direct_physical_f(:);
        h_source=exp(1i*k*cfg.s0_values_m(ss)).*h_plane;
        R=cfg.s0_values_m(ss)+cfg.L_values_m(ll);
        h_exact=cfg.green_amplitude.*exp(1i*k*R)./R;
        for ff=1:numel(f)
            q=q+1; row=local_axis_row(); row.s0_m=cfg.s0_values_m(ss);
            row.L_m=cfg.L_values_m(ll); row.range_m=R; row.frequency_hz=f(ff);
            row.H_plane_pe=h_plane(ff); row.H_source_pe=h_source(ff); row.H_exact=h_exact(ff);
            row.pe_tl_db=-20*log10(abs(h_source(ff)));
            row.exact_tl_db=-20*log10(abs(h_exact(ff)));
            row.tl_error_db=row.pe_tl_db-row.exact_tl_db;
            row.phase_error_rad=angle(h_source(ff)*conj(h_exact(ff)));
            axis_rows(q)=row;
        end
    end
end
axis_table=struct2table(axis_rows);

fine_f=cfg.f_ref_hz+cfg.group_delay_offsets_hz(:).';
p=local_params(cfg,fine_f,cfg.representative_L_m,cfg.representative_s0_m, ...
    0,cfg.axis_nx,cfg.axis_width_m,cfg.axis_stepz_lamb,'rx_only',cfg.axis_alpha_max_np_per_m);
fine=vertical_channel_model(p); fine_k=2*pi*fine.f_axis(:)/cfg.c0_mps;
fine_h=exp(1i*fine_k*cfg.representative_s0_m).*fine.H_direct_physical_f(:);
phase=unwrap(angle(fine_h)); slope=polyfit(fine.f_axis(:),phase,1);
group_delay_s=slope(1)/(2*pi);
exact_delay_s=(cfg.representative_s0_m+cfg.representative_L_m)/cfg.c0_mps;
else
    axis_table=struct2table(repmat(local_axis_row(),0,1));
    group_delay_s=NaN; exact_delay_s=NaN;
end

checks=table(["pe_as_full_field";"pe_as_trusted_field";"pe_as_phase"; ...
    "channel_invariant";"reflection_disabled";"dz_invariance"], ...
    [max(field_table.pe_as_rel_l2);max(field_table.pe_as_trusted_rel_l2); ...
    max(field_table.pe_as_phase_rms_rad);max(field_table.invariant_error); ...
    max(field_table.max_reflection);local_dz_invariance(field_table)], ...
    [cfg.pe_as_tolerance;cfg.pe_as_tolerance;cfg.pe_as_phase_tolerance_rad; ...
    cfg.invariant_tolerance;cfg.invariant_tolerance;cfg.dz_invariance_tolerance], ...
    'VariableNames',{'check_name','value','limit'});
checks.passed=checks.value<=checks.limit;

validation=struct('schema_version','1.0.0','config',cfg,'field_table',field_table, ...
    'field_data',{field_data},'axis_table',axis_table,'group_delay_s',group_delay_s, ...
    'exact_group_delay_s',exact_delay_s,'group_delay_error_s',group_delay_s-exact_delay_s, ...
    'checks',checks,'level1_passed',all(checks.passed));
validation.files=local_write_outputs(validation);
if ~validation.level1_passed
    error('validate_pe_as_freefield_vertical:Level1Failed', ...
        'PE-AS hard check failed; Bellhop comparison must not proceed.');
end
end

function cfg=local_defaults(root)
g(1)=struct('name','A_128_32_dz050','nx',128,'width_m',32,'stepz_lamb',0.5);
g(2)=struct('name','B_256_32_dz050','nx',256,'width_m',32,'stepz_lamb',0.5);
g(3)=struct('name','C_256_64_dz050','nx',256,'width_m',64,'stepz_lamb',0.5);
g(4)=struct('name','D_256_32_dz025','nx',256,'width_m',32,'stepz_lamb',0.25);
cfg=struct('output_dir',fullfile(root,'results','validation','pe_bellhop_freefield','pe_reference'), ...
    'c0_mps',1500,'f_ref_hz',4000,'frequencies_hz',[3000 4000 5000], ...
    'L_values_m',[20 40 70 100],'s0_values_m',[5 10 20], ...
    'representative_L_m',70,'representative_s0_m',10,'green_amplitude',1/(4*pi), ...
    'grid_cases',g,'axis_nx',256,'axis_width_m',64,'axis_stepz_lamb',0.5, ...
    'axis_alpha_max_np_per_m',0.15, ...
    'trusted_radius_m',2,'phase_floor_ratio',1e-6, ...
    'group_delay_offsets_hz',[-1 0 1],'pe_as_tolerance',1e-10, ...
    'pe_as_phase_tolerance_rad',1e-10,'invariant_tolerance',1e-12, ...
    'dz_invariance_tolerance',1e-10,'run_full_matrix',true);
end

function cfg=local_overrides(cfg,o)
names=fieldnames(o); for ii=1:numel(names), if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end, cfg.(names{ii})=o.(names{ii}); end
end

function p=local_params(cfg,f,L,s0,xrx,nx,width,step,save_mode,alpha)
p=struct('f0',f,'enable_wideband',false,'f_ref_hz',cfg.f_ref_hz,'c0',cfg.c0_mps, ...
    'z_max',L,'z_tx',L,'z_rx',0,'xw',width,'yw',width,'nx',nx,'ny',nx, ...
    'x_tx',0,'y_tx',0,'x_rx',xrx,'y_rx',0,'stepz_lamb',step, ...
    'source_mode','custom_field_fn','source_field_fn',@virtual_point_source_initial_field_vertical, ...
    'virtual_source_distance_m',s0,'source_green_amplitude',cfg.green_amplitude, ...
    'sponge_ratio',0.12,'alpha_max_np_per_m',alpha,'env_mode','uniform', ...
    'enable_surface_reflection',false,'enable_bubbles',false,'enforce_1_over_R',false, ...
    'show_figures',false,'save_mode',save_mode,'use_gpu',false,'channel_phase_reference','direct_dsp');
end

function r=local_field_row()
r=struct('case_name',"",'nx',NaN,'width_m',NaN,'stepz_lamb',NaN, ...
    'pe_as_rel_l2',NaN,'pe_as_max_abs',NaN,'pe_as_trusted_rel_l2',NaN, ...
    'pe_as_phase_rms_rad',NaN,'pe_exact_trusted_rel_l2',NaN, ...
    'pe_exact_phase_rms_rad',NaN,'pe_exact_tl_rms_db',NaN, ...
    'invariant_error',NaN,'max_reflection',NaN);
end

function r=local_axis_row()
r=struct('s0_m',NaN,'L_m',NaN,'range_m',NaN,'frequency_hz',NaN, ...
    'H_plane_pe',complex(NaN),'H_source_pe',complex(NaN),'H_exact',complex(NaN), ...
    'pe_tl_db',NaN,'exact_tl_db',NaN,'tl_error_db',NaN,'phase_error_rad',NaN);
end

function v=local_phase_rms(a,b,mask)
d=angle(a(mask).*conj(b(mask))); v=sqrt(mean(d.^2));
end

function v=local_tl_rms(a,b,mask)
d=-20*log10(max(abs(a(mask)),realmin))+20*log10(max(abs(b(mask)),realmin)); v=sqrt(mean(d.^2));
end

function value=local_dz_invariance(t)
value=0;
for ii=1:height(t)
    match=find(t.nx==t.nx(ii) & t.width_m==t.width_m(ii) & t.stepz_lamb~=t.stepz_lamb(ii));
    if ~isempty(match)
        value=max(value,max(abs(t.pe_as_rel_l2(ii)-t.pe_as_rel_l2(match))));
    end
end
end

function files=local_write_outputs(v)
out=v.config.output_dir; mat_file=fullfile(out,'pe_as_freefield_validation.mat');
field_csv=fullfile(out,'pe_as_field_convergence.csv'); axis_csv=fullfile(out,'pe_spherical_axis.csv');
checks_csv=fullfile(out,'pe_as_hard_checks.csv'); fig_file=fullfile(out,'pe_as_freefield.png');
writetable(v.field_table,field_csv); writetable(v.axis_table,axis_csv); writetable(v.checks,checks_csv);
fig=figure('Visible','off','Color','w','Position',[100 100 1100 750]); cleanup=onCleanup(@()close(fig));
d=v.field_data{min(2,numel(v.field_data))}; [~,iy]=min(abs(d.y_m));
subplot(2,2,1); plot(d.x_m,abs(d.psi_pe(iy,:)),'-',d.x_m,abs(d.psi_as(iy,:)),'--'); grid on; xlabel('x (m)'); ylabel('|Psi|'); legend('PE','AS'); title('Reduced field');
subplot(2,2,2); plot(d.x_m,unwrap(angle(d.psi_pe(iy,:))),'-',d.x_m,unwrap(angle(d.psi_as(iy,:))),'--'); grid on; xlabel('x (m)'); ylabel('phase (rad)'); title('Reduced phase');
subplot(2,2,3); imagesc(d.x_m,d.y_m,log10(max(abs(d.psi_pe-d.psi_as),1e-18))); axis image xy; colorbar; title('log10 |PE-AS|'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,2,4); semilogy(v.field_table.nx,v.field_table.pe_as_rel_l2,'o-'); grid on; xlabel('nx'); ylabel('relative L2'); title('PE-AS hard check');
exportgraphics(fig,fig_file,'Resolution',180); clear cleanup
schema_version=v.schema_version; validation=v; save(mat_file,'validation','schema_version','-v7.3');
files=struct('mat',mat_file,'field_csv',field_csv,'axis_csv',axis_csv,'checks_csv',checks_csv,'figure',fig_file);
end
