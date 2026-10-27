function receipt=validate_two_dimensional_s5a_weak_pe_sensitivity(frozen_source_path,output_dir,available_physical_bytes)
%VALIDATE_TWO_DIMENSIONAL_S5A_WEAK_PE_SENSITIVITY Same source, 6 kHz P22.
arguments
    frozen_source_path (1,:) char
    output_dir (1,:) char
    available_physical_bytes (1,1) double {mustBeFinite,mustBePositive}
end
assert(isfile(frozen_source_path) && ~isfolder(output_dir) && ~isfile(output_dir));
estimate=128*1024^2;
assert(estimate<=.7*available_physical_bytes,'RAM gate');
s=load(frozen_source_path,'src'); src=closeout2d_source_at_frequency(s.src,6000,1500);
N=984; W=192.1875;
x=(-W/2)+(0:N-1)'*W/N;
source=ifft(src.coeff(:));
[eta,~]=closeout_s2b_surface(x,.01);
cfg=struct('frequency_hz',6000,'c0_mps',1500,'xw_m',W,'nx',N, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
    'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
    'reflection_model','model0_normal','source_field_override',source, ...
    'eta_z_m',eta);
tic_run=tic;
base=run_pe_physical_surface_model0_validation(cfg);
half_cfg=cfg; half_cfg.step_m=.025;
half=run_pe_physical_surface_model0_validation(half_cfg);
fine_cfg=cfg; fine_cfg.nx=2*N;
xf=(-W/2)+(0:2*N-1)'*W/(2*N);
coeff_fine=complex(zeros(2*N,1));
coeff_fine(1:N/2)=2*src.coeff(1:N/2);
coeff_fine(end-N/2+1:end)=2*src.coeff(N/2+1:end);
fine_source=ifft(coeff_fine);
source_even_max_abs=max(abs(fine_source(1:2:end)-source));
assert(source_even_max_abs<1e-12,'Frozen-spectrum zero-padding changed source');
fine_cfg.source_field_override=fine_source;
[fine_eta,~]=closeout_s2b_surface(xf,.01);
fine_cfg.eta_z_m=fine_eta;
fine=run_pe_physical_surface_model0_validation(fine_cfg);
u_base=base.reflected_native(base.ix_rx);
u_half=half.reflected_native(half.ix_rx);
u_fine=fine.reflected_native(fine.ix_rx);
step_relative=abs(u_half-u_base)/abs(u_base);
grid_relative=abs(u_fine-u_base)/abs(u_base);
empirical_relative=max(step_relative,grid_relative);
metrics=table(6000,.01,step_relative,grid_relative,empirical_relative, ...
    source_even_max_abs,real(u_base),imag(u_base), ...
    'VariableNames',{'frequency_hz','amplitude_m','step_half_relative', ...
    'same_spectrum_double_grid_relative','empirical_relative', ...
    'source_even_max_abs','base_native_real','base_native_imag'});
receipt=struct('status','PE_SENSITIVITY_DIAGNOSTIC', ...
    'matlab_version',version,'frequency_hz',6000, ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimate,'pe_calls',3, ...
    'bie_boundary_solves',0,'elapsed_seconds',toc(tic_run));
mkdir(output_dir);
save(fullfile(output_dir,'pe_sensitivity.mat'),'metrics','receipt', ...
    'u_base','u_half','u_fine','source','fine_source');
writetable(metrics,fullfile(output_dir,'pe_sensitivity.csv'));
fid=fopen(fullfile(output_dir,'receipt.json'),'w');
assert(fid>0); guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt); disp(metrics);
end
