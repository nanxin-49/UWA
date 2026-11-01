function prediction = freeze_two_dimensional_s5_p22_halfstep_prediction(task_root)
%FREEZE_TWO_DIMENSIONAL_S5_P22_HALFSTEP_PREDICTION Blind P22 BIE holdout.
% Run only after all 4062.5 Hz BIE layers qualify, before any 4031.25 Hz
% BIE job. No acoustic boundary solve or fitted source/path phase is used.
arguments
    task_root (1,:) char
end
assert(isfolder(task_root) && isfolder(fullfile(task_root,'code')));
addpath(fullfile(task_root,'code'));
prediction_dir=fullfile(task_root,'prediction');
assert(~isfolder(prediction_dir) && ~isfile(prediction_dir), ...
    'Prediction already exists; refuse overwrite');
halfstep_dir=fullfile(task_root,'jobs','f4031p25');
if isfolder(halfstep_dir)
    ids=dir(halfstep_dir);
    assert(~any(ismember({ids.name},{'main96','self64','space14','window60','window70'})), ...
        'Half-step BIE job directory exists; prediction is no longer blind');
end
source_path=fullfile(task_root,'input','frozen_source.mat');
reference_path=fullfile(task_root,'input','fields_core.csv');
assert(isfile(source_path) && isfile(reference_path),'Frozen input missing');
source_sha=local_sha256(source_path);
reference_sha=local_sha256(reference_path);
assert(strcmp(source_sha, ...
    '8d60d3d9d6aa3aa16110ac33d81e0d193a61526aef2c0c2e9edd531cae531f55'), ...
    'Frozen source SHA mismatch');
assert(strcmp(reference_sha, ...
    '1b6b8656423ff39e2678639f52eb1be2748e75ad3e1786d6c05e2522f875ffd5'), ...
    '4 kHz reference CSV SHA mismatch');
s=load(source_path,'src'); frozen_src=s.src;
assert(frozen_src.n_original==984 && frozen_src.z_tx_m==100);
t=readtable(reference_path,'TextType','string');
ix=t.case_name=="weak_low_K" & t.method=="BIE" & ...
    t.receiver_id=="P22" & t.frequency_hz==4000;
assert(nnz(ix)==1,'4 kHz qualified P22 BIE reference row missing');
R4000=complex(t.u_native_real(ix),t.u_native_imag(ix));
assert(isfinite(R4000) && abs(R4000)>0);
ids={'main96','self64','space14','window60','window70'};
U=complex(zeros(5,1));
layer_sha=cell(5,1);
for j=1:5
    out=fullfile(task_root,'jobs','f4062p5',ids{j});
    result_path=fullfile(out,'result.mat');
    status_path=fullfile(out,'exit_code.txt');
    assert(isfile(result_path) && isfile(status_path),'Base-frequency layer incomplete');
    assert(str2double(strtrim(fileread(status_path)))==0,'Base-frequency MATLAB process failed');
    q=load(result_path,'boundary','receipt','pre');
    assert(strcmp(q.receipt.status,'COMPLETE_UNAUDITED') && ...
        strcmp(q.receipt.job_id,ids{j}) && q.receipt.frequency_hz==4062.5);
    assert(q.pre.frequency_hz==4062.5 && q.pre.surface_A_m==.01 && ...
        q.pre.surface_K_radpm==.1 && q.pre.surface_phi_rad==0 && ...
        q.pre.receiver_x_m==0 && q.pre.receiver_z_m==3 && ...
        q.pre.tau_ref_s==96.875/1500);
    assert(q.boundary.finite && q.boundary.frequency_hz==4062.5 && ...
        q.boundary.receiver_x_m==0 && q.boundary.receiver_z_m==3 && ...
        q.boundary.node_count==q.pre.node_estimate);
    U(j)=q.boundary.receiver_field;
    assert(isfinite(U(j)) && abs(U(j))>0);
    assert(abs(U(j)-complex(q.receipt.receiver_native_real, ...
        q.receipt.receiver_native_imag))<1e-14);
    layer_sha{j}=local_sha256(result_path);
end
reference_sensitivity=max(abs(U(2:end)-U(1)))/abs(U(1));
assert(reference_sensitivity<=1e-4, ...
    '4062.5 Hz BIE empirical reference gate failed');
f=[4000 4062.5 4031.25];
F=complex(zeros(1,3));
for j=1:3
    src=closeout2d_source_at_frequency(frozen_src,f(j),1500);
    F(j)=-fk_single_receiver_incident(src,0,-3);
    assert(isfinite(F(j)) && abs(F(j))>0,'Flat reference invalid');
end
Q4000=R4000/F(1);
Q4062=U(1)/F(2);
Qpred=(Q4000+Q4062)/2;
Rpred=F(3)*Qpred;
assert(isfinite(Rpred));
prediction=struct('status','FROZEN_BEFORE_HALFSTEP_BIE', ...
    'frequency_hz',4031.25,'base_frequencies_hz',[4000,4062.5], ...
    'native_real',real(Rpred),'native_imag',imag(Rpred), ...
    'R4000_real',real(R4000),'R4000_imag',imag(R4000), ...
    'R4062p5_real',real(U(1)),'R4062p5_imag',imag(U(1)), ...
    'F_real',real(F),'F_imag',imag(F), ...
    'Q4000_real',real(Q4000),'Q4000_imag',imag(Q4000), ...
    'Q4062p5_real',real(Q4062),'Q4062p5_imag',imag(Q4062), ...
    'reference_sensitivity_relative',reference_sensitivity, ...
    'source_sha256',source_sha,'reference_csv_sha256',reference_sha, ...
    'base_layer_sha256',{layer_sha}, ...
    'formula','Rpred=F4031.25*(R4000/F4000+R4062.5/F4062.5)/2', ...
    'phase','native exp(-i omega t); no amplitude or phase fit', ...
    'tau_ref_s',96.875/1500, ...
    'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
mkdir(prediction_dir);
json_path=fullfile(prediction_dir,'prediction.json');
mat_path=fullfile(prediction_dir,'prediction.mat');
save(mat_path,'prediction','Rpred','R4000','U','F','Q4000','Q4062','Qpred');
fid=fopen(json_path,'w'); assert(fid>0);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(prediction));
clear guard
fid=fopen(fullfile(prediction_dir,'prediction_sha256.txt'),'w'); assert(fid>0);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',local_sha256(json_path));
disp(prediction);
end

function h=local_sha256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    md.update(fread(fid,1048576,'*uint8'));
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
