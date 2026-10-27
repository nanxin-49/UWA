function receipt=prepare_two_dimensional_s5a_8k_weak_local(frozen_source_path,output_dir,available_physical_bytes)
%PREPARE_TWO_DIMENSIONAL_S5A_8K_WEAK_LOCAL Bounded A-stage P22 model controls.
% No BIE reference or S5 broadband claim follows from this diagnostic.
arguments
    frozen_source_path (1,:) char
    output_dir (1,:) char
    available_physical_bytes (1,1) double {mustBeFinite,mustBePositive}
end
assert(isfile(frozen_source_path) && ~isfolder(output_dir) && ~isfile(output_dir));
estimate=256*1024^2;
assert(estimate<=.7*available_physical_bytes,'RAM gate');
s=load(frozen_source_path,'src'); src=closeout2d_source_at_frequency(s.src,8000,1500);
N=984; W=192.1875; A=.01; f=8000; c0=1500;
x=(-W/2)+(0:N-1)'*W/N;
source_tx=ifft(src.coeff(:));
[eta,~]=closeout_s2b_surface(x,A);
cfg=struct('frequency_hz',f,'c0_mps',c0,'xw_m',W,'nx',N, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
    'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
    'reflection_model','model0_normal','source_field_override',source_tx, ...
    'eta_z_m',eta);
tic_run=tic;
pe=run_pe_physical_surface_model0_validation(cfg);
half_cfg=cfg; half_cfg.step_m=.025;
pe_half=run_pe_physical_surface_model0_validation(half_cfg);
fine_cfg=cfg; fine_cfg.nx=2*N;
xf=(-W/2)+(0:2*N-1)'*W/(2*N);
coeff_fine=complex(zeros(2*N,1));
coeff_fine(1:N/2)=2*src.coeff(1:N/2);
coeff_fine(end-N/2+1:end)=2*src.coeff(N/2+1:end);
fine_source=ifft(coeff_fine);
assert(max(abs(fine_source(1:2:end)-source_tx))<1e-12);
fine_cfg.source_field_override=fine_source;
[fine_eta,~]=closeout_s2b_surface(xf,A);
fine_cfg.eta_z_m=fine_eta;
pe_fine=run_pe_physical_surface_model0_validation(fine_cfg);
u_direct=fk_single_receiver_incident(src,0,3);
u_flat=-fk_single_receiver_incident(src,0,-3);
u_pe=pe.reflected_native(pe.ix_rx);
u_half=pe_half.reflected_native(pe_half.ix_rx);
u_fine=pe_fine.reflected_native(pe_fine.ix_rx);
direct_pe=pe.direct_receiver*exp(1i*2*pi*f/c0*97);
direct_relative=abs(direct_pe-u_direct)/abs(u_direct);
assert(direct_relative<1e-10,'Direct source field control');
fk400=closeout_s2b_fullfk_analytic_exterior(src,A,0,3,400,u_flat);
fk800=closeout_s2b_fullfk_analytic_exterior(src,A,0,3,800,u_flat);
step_relative=abs(u_half-u_pe)/abs(u_pe);
grid_relative=abs(u_fine-u_pe)/abs(u_pe);
fk_quadrature_relative=abs(fk400.receiver_field-fk800.receiver_field)/ ...
    abs(fk800.receiver_field);
phase=exp(1i*2*pi*f*97/c0);
H_direct=conj(u_direct)*phase;
H_PE_reflect=conj(u_pe)*phase;
H_FullFK_reflect=conj(fk800.receiver_field)*phase;
metrics=table(f,A,direct_relative,step_relative,grid_relative, ...
    fk_quadrature_relative,real(H_direct),imag(H_direct), ...
    real(H_PE_reflect),imag(H_PE_reflect), ...
    real(H_FullFK_reflect),imag(H_FullFK_reflect), ...
    'VariableNames',{'frequency_hz','amplitude_m','direct_relative', ...
    'pe_halfstep_relative','pe_same_spectrum_grid_relative', ...
    'fk_400_800_relative','H_direct_real','H_direct_imag', ...
    'H_PE_reflect_real','H_PE_reflect_imag', ...
    'H_FullFK_reflect_real','H_FullFK_reflect_imag'});
receipt=struct('status','PRE_REFERENCE_DIAGNOSTIC_ONLY', ...
    'matlab_version',version,'frequency_hz',f,'amplitude_m',A, ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimate,'pe_calls',3, ...
    'fullfk_evaluations',2,'bie_boundary_solves',0, ...
    'elapsed_seconds',toc(tic_run));
mkdir(output_dir);
save(fullfile(output_dir,'model_controls.mat'),'metrics','receipt','src', ...
    'u_direct','u_flat','u_pe','u_half','u_fine','fk400','fk800');
writetable(metrics,fullfile(output_dir,'model_controls.csv'));
fid=fopen(fullfile(output_dir,'receipt.json'),'w');
assert(fid>0); guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt); disp(metrics(:,1:6));
end
