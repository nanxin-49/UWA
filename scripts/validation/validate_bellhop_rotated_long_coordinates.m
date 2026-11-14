function receipt = validate_bellhop_rotated_long_coordinates(overlay_exe,output_dir)
%VALIDATE_BELLHOP_ROTATED_LONG_COORDINATES Synthetic precision stress only.
% Six ranges over 5cm at 70km; a uniform field with 301 beams is inexpensive.
% This is not a 70km Munk propagation qualification or a production scenario.
arguments
    overlay_exe (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM gate.');
cfg=make_bellhop_rotated_munk_config('up'); cfg.uniform_control=true; cfg.point3d=true;
cfg.receiver_range_m=70000+(0:.01:.05); cfg.receiver_q_m=[-1 0 1];
cfg.range_limit_m=70000.1; cfg.beam_count=301; cfg.angle_limits_deg=[-.001 .001]; cfg.step_m=100;
stem=fullfile(output_dir,'case'); meta=write_bellhop_rotated_horizontal_q(stem,cfg,1);
process=run_bellhop_rotated_local_case(overlay_exe,stem,15);
field=read_bellhop_rotated_field(stem);
distance=hypot(cfg.receiver_q_m(:),cfg.receiver_range_m);
exact=exp(1i*2*pi*10000*distance/1500)./distance;
p=-conj(field.pressure);
complex_error=norm(p(:)-exact(:))/norm(exact(:));
phase_error=max(abs(angle(p.*conj(exact))),[],'all');
phase_increment=angle(p(2,2:end).*conj(p(2,1:end-1)));
expected_increment=2*pi*10000*diff(cfg.receiver_range_m)/1500;
increment_error=max(abs(angle(exp(1i*(phase_increment-expected_increment)))));
checks=table(["long_range_uniform_Green";"long_range_phase_rad";"one_cm_phase_increment_rad";"all_actual_ranges_distinct"], ...
    [complex_error<=1e-6;phase_error<=1e-6;increment_error<=1e-6;all(diff(field.receiver_range_m)>0)], ...
    [complex_error;phase_error;increment_error;double(any(diff(field.receiver_range_m)<=0))], ...
    [1e-6;1e-6;1e-6;0],'VariableNames',{'check','pass','value','limit'});
receipt=struct('status','COMPLETE','failed_checks',nnz(~checks.pass),'elapsed_s',toc(t), ...
    'matlab_release',version,'result_dir',output_dir, ...
    'scope','Synthetic float64 coordinate/phase stress, not 70km Munk field', ...
    'metadata_max_error_m',max(abs(field.receiver_range_m-field.receiver_range_metadata_m)), ...
    'overlay_exe_sha256',black_file_hash(overlay_exe));
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'cfg','meta','process','field','exact','p','checks','receipt');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
end
