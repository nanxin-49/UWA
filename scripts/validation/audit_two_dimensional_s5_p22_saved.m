function audit = audit_two_dimensional_s5_p22_saved(output_dir)
%AUDIT_TWO_DIMENSIONAL_S5_P22_SAVED Recheck scalar saved pilot evidence only.
% No propagation, field synthesis, remote access or function-handle loading.
% SHA manifests are checked separately by the baseline repository receipt.
arguments
    output_dir (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(java.io.File(output_dir).isAbsolute());
assert(~isfolder(output_dir) && ~isfile(output_dir));
mkdir(output_dir);
base=fullfile(root,'results','validation','two_dimensional_s5_p22_dualfreq', ...
    '20261031_01');
mirror=fullfile(base,'audit_mirror');
historical=jsondecode(fileread(fullfile(base,'audit_dualfreq_complete_01','audit.json')));
prediction=jsondecode(fileread(fullfile(mirror,'prediction','prediction.json')));
frozen=load(fullfile(mirror,'input','frozen_source.mat'),'src');
assert(frozen.src.n_original==984 && frozen.src.z_tx_m==100);
layers={'main96','self64','space14','window60','window70'};
labels={'f4062p5','f4031p25'}; frequencies=[4062.5,4031.25];
tau_ref=96.875/1500;
sensitivity=zeros(1,2); fk_error=sensitivity; pe_error=sensitivity;
fk_quadrature=sensitivity; references=complex(sensitivity);
native_fields=complex(zeros(5,2));
for k=1:2
    f=frequencies(k); phase=exp(1i*2*pi*f*tau_ref);
    jobs=fullfile(mirror,'jobs',labels{k});
    model=fullfile(jobs,'model','result.mat');
    assert(scalar(model,'/cfg/frequency_hz')==f);
    assert(scalar(model,'/cfg/step_m')==.05 && scalar(model,'/cfg/nx')==984);
    assert(scalar(model,'/tau_ref')==tau_ref);
    assert_source(model,frozen.src);
    direct=scalar(model,'/H_direct_native');
    close_value(scalar(model,'/H_pe_total_native'),direct+scalar(model,'/H_pe_reflect_native'));
    close_value(scalar(model,'/H_fk_total_native'),direct+scalar(model,'/H_fk_reflect_native'));
    close_value(scalar(model,'/H_pe_dsp'),conj(scalar(model,'/H_pe_total_native'))*phase);
    close_value(scalar(model,'/H_fk_dsp'),conj(scalar(model,'/H_fk_total_native'))*phase);
    fk=scalar(model,'/fk800/receiver_field');
    close_value(fk,scalar(model,'/H_fk_reflect_native'));
    fk_quadrature(k)=abs(scalar(model,'/fk400/receiver_field')-fk)/abs(fk);
    mr=job_receipt(fullfile(jobs,'model'),f);
    close_value(fk_quadrature(k),mr.fk_400_800_relative);
    for j=1:5
        job=fullfile(jobs,layers{j}); path=fullfile(job,'result.mat');
        r=job_receipt(job,f);
        assert(strcmp(r.job_id,layers{j}));
        assert(scalar(path,'/pre/frequency_hz')==f);
        assert(scalar(path,'/boundary/frequency_hz')==f);
        assert(scalar(path,'/pre/tau_ref_s')==tau_ref);
        assert(scalar(path,'/tau_ref')==tau_ref);
        assert(scalar(path,'/pre/receiver_x_m')==0 && scalar(path,'/pre/receiver_z_m')==3);
        assert(scalar(path,'/pre/surface_A_m')==.01);
        assert(scalar(path,'/pre/surface_K_radpm')==.1 && scalar(path,'/pre/surface_phi_rad')==0);
        assert(scalar(path,'/boundary/finite')==1);
        assert(scalar(path,'/boundary/node_count')==scalar(path,'/pre/node_estimate'));
        assert_source(path,frozen.src);
        value=scalar(path,'/H_reflect_native'); native_fields(j,k)=value;
        close_value(value,scalar(path,'/boundary/receiver_field'),1e-14);
        close_value(value,complex(r.receiver_native_real,r.receiver_native_imag),1e-14);
        close_value(scalar(path,'/H_direct_native'),direct);
        close_value(scalar(path,'/H_total_native'),direct+value);
        close_value(scalar(path,'/H_reflect_dsp'),conj(value)*phase);
        close_value(scalar(path,'/H_total_dsp'),conj(direct+value)*phase);
        close_value(scalar(path,'/boundary/linear_residual'),r.linear_residual,1e-14);
        close_value(scalar(path,'/boundary/offgrid_boundary_residual'),r.offgrid_boundary_residual,1e-14);
    end
    references(k)=native_fields(1,k);
    sensitivity(k)=max(abs(native_fields(2:end,k)-references(k)))/abs(references(k));
    assert(sensitivity(k)<=1e-4);
    pe_error(k)=abs(scalar(model,'/H_pe_reflect_native')-references(k))/abs(references(k));
    fk_error(k)=abs(fk-references(k))/abs(references(k));
end
F=complex(prediction.F_real,prediction.F_imag);
pred=F(3)*(complex(prediction.R4000_real,prediction.R4000_imag)/F(1)+ ...
    complex(prediction.R4062p5_real,prediction.R4062p5_imag)/F(2))/2;
saved=complex(prediction.native_real,prediction.native_imag);
close_value(pred,saved,1e-14);
assert(strcmp(prediction.status,'FROZEN_BEFORE_HALFSTEP_BIE'));
assert(prediction.frequency_hz==4031.25 && prediction.tau_ref_s==tau_ref);
prediction_error=abs(pred-references(2))/abs(references(2));
close_value(prediction_error,historical.prediction_complex_relative_error,1e-12);
close_value(sensitivity(1),historical.bie_reference_4062p5_relative);
close_value(sensitivity(2),historical.bie_reference_4031p25_relative);
resolved=fk_error>5*(sensitivity+fk_quadrature);
assert(~any(resolved));
audit=struct('status','SAVED_SCALAR_REFERENCE_RECHECK_PASS_SCOPED', ...
    'frequencies_hz',frequencies,'tau_ref_s',tau_ref, ...
    'bie_reference_relative',sensitivity,'bie_limit',1e-4, ...
    'pe_bie_relative',pe_error,'fk_bie_relative',fk_error, ...
    'fk_quadrature_relative',fk_quadrature,'fk_resolved_over_five',resolved, ...
    'prediction_complex_relative_error',prediction_error, ...
    'failed_checks',0,'new_solver_calls',0, ...
    'limitations',{{'Scalar fields, source arrays, configuration, exit records and algebra only.', ...
    'Historical source_at field synthesis and every saved boundary array not replayed.', ...
    'Prediction chronology remains historical log evidence; no independent clock attestation.', ...
    'PE step/grid sensitivity absent; no full-band CIR or communication qualification.'}});
save(fullfile(output_dir,'saved_fields.mat'),'native_fields','references','audit');
fid=fopen(fullfile(output_dir,'audit.json'),'w'); assert(fid>0);
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
disp(audit);
end

function value=scalar(path,key)
value=h5read(path,key);
if isstruct(value), value=complex(value.real,value.imag); end
assert(isscalar(value) && isfinite(value));
end

function close_value(a,b,tolerance)
if nargin<3, tolerance=1e-12; end
assert(abs(a-b)<=tolerance*max(1,abs(b)));
end

function assert_source(path,source)
kx=h5read(path,'/src/kx'); c=h5read(path,'/src/coeff');
if isstruct(c), c=complex(c.real,c.imag); end
assert(isequal(kx(:),source.kx(:)) && isequal(c(:),source.coeff(:)));
end

function receipt=job_receipt(job,frequency)
assert(str2double(strtrim(fileread(fullfile(job,'exit_code.txt'))))==0);
seconds=str2double(strtrim(fileread(fullfile(job,'process_elapsed_seconds.txt'))));
assert(seconds>0 && seconds<900);
assert(isfile(fullfile(job,'stdout.log')) && isfile(fullfile(job,'stderr.log')));
receipt=jsondecode(fileread(fullfile(job,'receipt.json')));
assert(receipt.frequency_hz==frequency && strcmp(receipt.status,'COMPLETE_UNAUDITED'));
end
