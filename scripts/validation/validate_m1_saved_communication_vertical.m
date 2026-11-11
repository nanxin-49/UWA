function receipt=validate_m1_saved_communication_vertical(run_dir)
%VALIDATE_M1_SAVED_COMMUNICATION_VERTICAL Replay existing M0 wideband input.
% Interface regression only: no new frequency solve, CIR study or comm stage.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
job=fullfile(run_dir,'saved_comm_control'); assert(~isfolder(job)); mkdir(job);
t=tic; old=maxNumCompThreads(1); c=onCleanup(@()maxNumCompThreads(old)); %#ok<NASGU>
diary(fullfile(job,'stdout.log')); d=onCleanup(@()diary('off')); %#ok<NASGU>
input=fullfile(root,'results','validation','m0_phase_contract','20261110_m0_native_dsp_v2_01','E_cir_comm_boundary_fixed','payload.mat');
s=load(input,'payload'); original=s.payload.data.production; saved=s.payload.data.production_comm;
set=struct('H_fm',original.H_f,'f_axis_hz',original.f_axis(:),'phase_reference_meta',original.phase_reference_meta);
output=evaluate_mpsk_channel_ensemble_vertical(set,struct('M',4,'n_sym',128, ...
    'tap_fft_size',64,'symbol_rate_hz',1000,'f_ref_hz',6000,'EbN0_dB_list',[], ...
    'bootstrap_count',10,'mmse_reg_eps',1e-12));
assert(output.noiseless_BER==0 && output.noiseless_SER==0);
assert(output.noiseless_BER==saved.noiseless_BER);
assert(strcmp(output.receive_window_mode,'full_linear'));
old_pilot=load(fullfile(root,'results','validation','m1_convergence','20261111_m1_execution_02','P4','payload.mat'),'payload');
new_pilot=load(fullfile(run_dir,'P4','payload.mat'),'payload');
names={'native_direct','native_reflect','native_total','physical_direct','physical_reflect','physical_total','dsp_direct','dsp_reflect','dsp_total'};
pilot_repeat_max_difference=0;
assert(isequal(old_pilot.payload.coords_m,new_pilot.payload.coords_m));
for j=1:numel(names)
    pilot_repeat_max_difference=max(pilot_repeat_max_difference,max(abs(old_pilot.payload.(names{j})-new_pilot.payload.(names{j}))));
end
assert(pilot_repeat_max_difference==0,'Repeated pilot changed numerical fields.');
clear old_pilot new_pilot
save(fullfile(job,'payload.mat'),'output','original','saved','pilot_repeat_max_difference');
receipt=struct('status','CONTROL_PASS','matlab_pid',feature('getpid'),'matlab_version',version, ...
    'elapsed_s',toc(t),'noiseless_BER',output.noiseless_BER,'noiseless_SER',output.noiseless_SER, ...
    'failed_checks',0,'input_path',input,'input_sha256',local_sha(input), ...
    'pilot_repeat_max_difference',pilot_repeat_max_difference, ...
    'code_sha256',local_sha([mfilename('fullpath') '.m']), ...
    'scope','Existing M0 saved H replay only; no new propagation, no formal communication qualification.');
b=jsondecode(fileread(fullfile(run_dir,'budget.json'))); b.used_s=b.used_s+receipt.elapsed_s;
assert(b.used_s<=14400);
fid=fopen(fullfile(run_dir,'budget.json'),'w'); fprintf(fid,'%s',jsonencode(b,PrettyPrint=true)); fclose(fid);
clear d
fid=fopen(fullfile(job,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('M1_SAVED_COMM_CONTROL %s\n',jsonencode(receipt));
end
function value=local_sha(file)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(file,'rb'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid), bytes=fread(fid,2^20,'*uint8'); if ~isempty(bytes), md.update(typecast(bytes,'int8')); end, end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
