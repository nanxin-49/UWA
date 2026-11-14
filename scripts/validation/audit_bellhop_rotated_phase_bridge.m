function receipt = audit_bellhop_rotated_phase_bridge(input_dir,output_dir,influence_source)
%AUDIT_BELLHOP_ROTATED_PHASE_BRIDGE Audit the fixed source phasor convention.
% Native geometric point output has ScalePressure const=-1; exp(-i*omega*tau).
% FFP exp(-i*omega*t) with positive Green requires -conj(p_Bellhop).
% The sign/conjugation is prescribed from source and uniform Green, not fitted.
arguments
    input_dir (1,:) char
    output_dir (1,:) char
    influence_source (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
source=fileread(influence_source);
assert(contains(source,"const = -1.0") && contains(source,"factor = const / SQRT( ABS( r( ir ) ) )"));
data=load(fullfile(input_dir,'result.mat'),'fields','receipt');
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
pressure_ffp_convention=cell(1,2); diagnostic_ratio=cell(1,2);
for j=1:2
    field=data.fields{j}; r=field.cfg.receiver_range_m;
    pressure_ffp_convention{j}=-conj(field.official.pressure);
    expected=field.native_transport_amplitude.*exp(1i*2*pi*field.cfg.frequency_hz*field.reference_time_s);
    phase_error=max(abs(angle(pressure_ffp_convention{j}.*conj(expected))));
    add(sprintf('field%d_fixed_bridge_phase_rad',j),phase_error,.05);
    add(sprintf('field%d_native_transport_amplitude',j), ...
        max(abs(abs(pressure_ffp_convention{j})./abs(expected)-1)),2e-3);
    if j==1
        exact=exp(1i*2*pi*field.cfg.frequency_hz*r/field.cfg.uniform_c_mps)./r;
        add('uniform_positive_Green_complex_relative',norm(pressure_ffp_convention{j}-exact)/norm(exact),2e-3);
    end
    diagnostic_ratio{j}=field.physical_to_native_amplitude_ratio;
end
receipt=struct('status','COMPLETE','elapsed_s',toc(t),'failed_checks',nnz(~checks.pass), ...
    'fixed_bridge','pressure_ffp_convention = -conj(pressure_Bellhop_geometric_point)', ...
    'source_sha256',black_file_hash(influence_source), ...
    'input_sha256',black_file_hash(fullfile(input_dir,'result.mat')), ...
    'scope','Source phasor bridge only. Munk native R amplitude still lacks rotated 3D spreading qualification.', ...
    'matlab_release',version,'result_dir',output_dir);
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'receipt','checks','pressure_ffp_convention','diagnostic_ratio');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end