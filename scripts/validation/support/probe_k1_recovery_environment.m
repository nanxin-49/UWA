function receipt=probe_k1_recovery_environment(root,mode)
%PROBE_K1_RECOVERY_ENVIRONMENT Nonacoustic A1/A2, no project physics on path.
arguments
    root (1,:) char
    mode (1,:) char {mustBeMember(mode,{'A1','A2'})}
end
assert(startsWith(root,'/home/misc/cwx/k1_recovery_diag_'));
folder=fullfile(root,mode); assert(isfolder(folder));
diary(fullfile(folder,'diary.log')); guard=onCleanup(@()diary('off'));
t=tic; event(root,mode,'BEGIN');
receipt=struct('mode',mode,'MATLAB_version',version,'PID',feature('getpid'), ...
    'single_computational_thread',maxNumCompThreads==1,'new_acoustic_solves',0,'failed_checks',0);
assert(strcmp(version('-release'),'2025b') && maxNumCompThreads==1);
if strcmp(mode,'A1')
    event(root,mode,'VERSION_PATH_LICENSE');
    receipt.installed_products=ver; receipt.basic_license_test=logical(license('test','MATLAB'));
    assert(receipt.basic_license_test);
    receipt.fft_resolution=which('fft'); receipt.readtable_resolution=which('readtable');
    receipt.BLAS=version('-blas');
    token=fullfile(folder,'own_log_roundtrip.txt'); fid=fopen(token,'w'); assert(fid>=0);
    fprintf(fid,'K1_A1_NONACOUSTIC_ONLY\n'); fclose(fid);
    assert(strcmp(strtrim(fileread(token)),'K1_A1_NONACOUSTIC_ONLY'));
else
    event(root,mode,'READ_FROZEN_N40_INPUT_NO_SURFACE_GENERATION');
    T=readtable(fullfile(root,'input','coefficients.csv'));
    cfg=jsondecode(fileread(fullfile(root,'input','frozen_config.json')));
    assert(height(T)==40 && cfg.modes==40 && cfg.U_PM_mps==6 && cfg.seed==260001);
    assert(all(isfinite(T{:,:}),'all'));
    receipt.N40_rows=height(T); receipt.input_only_no_surface_evaluation=true;
    event(root,mode,'SYNTHETIC_FFT');
    rng(260001,'twister'); x=randn(32768,1)+1i*randn(32768,1);
    y=ifft(fft(x)); receipt.FFT_roundtrip_relative=norm(x-y)/norm(x);
    assert(receipt.FFT_roundtrip_relative<1e-12);
    event(root,mode,'SYNTHETIC_LINEAR_ALGEBRA');
    A=randn(256); A=A'*A+eye(256); b=randn(256,1); z=A\b;
    receipt.linear_residual_relative=norm(A*z-b)/norm(b);
    assert(receipt.linear_residual_relative<1e-10);
    event(root,mode,'SYNTHETIC_SAVE_LOAD');
    save(fullfile(folder,'synthetic_roundtrip.mat'),'x','z','-v7');
    q=load(fullfile(folder,'synthetic_roundtrip.mat')); assert(isequal(q.x,x) && isequal(q.z,z));
    receipt.synthetic_only=true;
end
event(root,mode,'END'); receipt.status='NONACOUSTIC_DIAGNOSTIC_PASS'; receipt.elapsed_s=toc(t);
fid=fopen(fullfile(folder,'receipt.json'),'w'); assert(fid>=0); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('K1_%s_NONACOUSTIC_PASS elapsed=%.6f\n',mode,receipt.elapsed_s);
end
function event(root,mode,label)
fid=fopen(fullfile(root,mode,'events.log'),'a'); assert(fid>=0);
fprintf(fid,'%s %s pid=%d\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS')),label,feature('getpid')); fclose(fid);
fprintf('K1_DIAGNOSTIC %s %s\n',mode,label);
end
