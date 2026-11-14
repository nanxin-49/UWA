function receipt = validate_bellhop_rotated_disabled_mode(reference_dir,overlay_exe,output_dir)
%VALIDATE_BELLHOP_ROTATED_DISABLED_MODE Reconstruct one frozen direct-Q case.
arguments
    reference_dir (1,:) char
    overlay_exe (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM gate.');
saved=load(fullfile(reference_dir,'result.mat'),'fields');
cfg=saved.fields{2}.cfg;
stem=fullfile(output_dir,'case'); write_bellhop_rotated_munk_q(stem,cfg);
source=fullfile(reference_dir,'field2_binary1','case');
assert(strcmp(black_file_hash([stem '.env']),black_file_hash([source '.env'])),'Disabled ENV bytes changed.');
assert(strcmp(black_file_hash([stem '.ssp']),black_file_hash([source '.ssp'])),'Disabled Q bytes changed.');
fid=fopen([stem '.iwcov'],'w'); assert(fid>=0); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
process=run_bellhop_rotated_local_case(overlay_exe,stem,20);
actual=read_bellhop_shd_unfolded_vertical([stem '.shd']);
old=saved.fields{2}.official.pressure;
relative=norm(actual.pressure(:)-old(:))/norm(old(:));
checks=table("disabled_direct_complex_relative",relative<=1e-4,relative,1e-4, ...
    'VariableNames',{'check','pass','value','limit'});
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status','COMPLETE','failed_checks',nnz(~checks.pass),'elapsed_s',toc(t), ...
    'overlay_exe_sha256',black_file_hash(overlay_exe),'result_dir',output_dir, ...
    'scope','Current writer default ENV/SSP byte-identical; direct-Q complex pressure vs official frozen case', ...
    'matlab_release',version);
save(fullfile(output_dir,'result.mat'),'receipt','checks','actual','process','cfg');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
end