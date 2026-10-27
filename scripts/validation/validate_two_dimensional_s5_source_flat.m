function receipt = validate_two_dimensional_s5_source_flat(frozen_source_path,output_dir,available_physical_bytes)
%VALIDATE_TWO_DIMENSIONAL_S5_SOURCE_FLAT Bounded three-frequency P22 flat control.
% No BIE boundary solve or rough-surface solve is performed.
arguments
    frozen_source_path (1,:) char
    output_dir (1,:) char
    available_physical_bytes (1,1) double {mustBeFinite,mustBePositive}
end
assert(isfile(frozen_source_path) && ~isfolder(output_dir) && ~isfile(output_dir));
estimated_increment_bytes=128*1024^2;
assert(estimated_increment_bytes<=0.7*available_physical_bytes,'RAM gate');
s=load(frozen_source_path,'src'); frozen=s.src;
N=frozen.n_original; assert(N==984);
W=192.1875; x0=-W/2;
x=x0+(0:N-1)'*W/N;
assert(abs(frozen.x0_m-x0)<1e-12 && frozen.z_tx_m==100);
source_tx=ifft(frozen.coeff(:));
assert(max(abs(source_tx-exp(-0.5*(x/0.3).^2)))<1e-12);
frequencies_hz=[4000;6000;8000];
flat_rows=cell(numel(frequencies_hz),1);
tic_total=tic;
for ii=1:numel(frequencies_hz)
    f=frequencies_hz(ii); k=2*pi*f/1500;
    src=closeout2d_source_at_frequency(frozen,f,1500);
    if ii==1
        assert(norm(src.kz-frozen.kz)/norm(frozen.kz)<1e-14,'4 kHz anchor changed');
    end
    tx=fk_single_receiver_incident(src,x,100+zeros(N,1));
    source_plane_error=norm(tx-source_tx)/norm(source_tx);
    direct_native=fk_single_receiver_incident(src,0,3);
    reflected_native=-fk_single_receiver_incident(src,0,-3);
    cfg=struct('frequency_hz',f,'c0_mps',1500,'xw_m',W,'nx',N, ...
        'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3, ...
        'surface_reflect_coeff',-1,'step_m',0.05,'x_rx_m',0, ...
        'reflection_model','model0_normal','source_field_override',source_tx, ...
        'eta_z_m',zeros(N,1));
    tic_pe=tic; pe=run_pe_physical_surface_model0_validation(cfg);
    pe_seconds=toc(tic_pe);
    pe_direct_native=pe.direct_receiver*exp(1i*k*97);
    pe_reflected_native=pe.reflected_native(pe.ix_rx);
    direct_error=abs(pe_direct_native-direct_native)/abs(direct_native);
    reflected_error=abs(pe_reflected_native-reflected_native)/abs(reflected_native);
    tau_d=97/1500;
    H_direct=conj(direct_native)*exp(1i*2*pi*f*tau_d);
    H_reflect=conj(reflected_native)*exp(1i*2*pi*f*tau_d);
    H_total=H_direct+H_reflect;
    pe_H_total=conj(pe_direct_native+pe_reflected_native)*exp(1i*2*pi*f*tau_d);
    total_error=abs(pe_H_total-H_total)/abs(H_total);
    assert(source_plane_error<1e-10 && direct_error<1e-10 && ...
        reflected_error<1e-10 && total_error<1e-10,'Flat source/field gate');
    flat_rows{ii}=table(f,source_plane_error,direct_error,reflected_error, ...
        total_error,real(H_direct),imag(H_direct),real(H_reflect), ...
        imag(H_reflect),real(H_total),imag(H_total),pe_seconds, ...
        'VariableNames',{'frequency_hz','tx_source_relative', ...
        'pe_direct_relative','pe_reflected_relative','pe_total_relative', ...
        'H_direct_real','H_direct_imag','H_reflect_real','H_reflect_imag', ...
        'H_total_real','H_total_imag','pe_wall_seconds'});
end
metrics=vertcat(flat_rows{:});
receipt=struct('status','THREE_FREQUENCY_FLAT_P22_PASS_ONLY', ...
    'matlab_version',version,'frequencies_hz',frequencies_hz, ...
    'source_kind','frozen 984-term Gaussian Tx spatial field', ...
    'field_unit','unit-source response, not Pa', ...
    'phase_convention','native exp(-i omega t)', ...
    'common_dsp_origin_s',97/1500,'estimated_increment_bytes', ...
    estimated_increment_bytes,'available_physical_bytes',available_physical_bytes, ...
    'pe_calls',numel(frequencies_hz),'bie_boundary_solves',0, ...
    'fullfk_rough_integrals',0,'elapsed_seconds',toc(tic_total));
mkdir(output_dir);
save(fullfile(output_dir,'flat_controls.mat'),'metrics','receipt','frequencies_hz', ...
    'source_tx','frozen_source_path');
writetable(metrics,fullfile(output_dir,'flat_metrics.csv'));
fid=fopen(fullfile(output_dir,'receipt.json'),'w');
assert(fid>0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt); disp(metrics(:,1:5));
end
