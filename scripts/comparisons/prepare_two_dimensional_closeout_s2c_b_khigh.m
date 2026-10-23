function p = prepare_two_dimensional_closeout_s2c_b_khigh(run_id,available_physical_bytes)
%PREPARE_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_KHIGH P0/P1, no acoustic evaluation.
arguments
    run_id (1,:) char
    available_physical_bytes (1,1) double {mustBePositive}
end
assert(~isempty(regexp(run_id,'^20261023_s2c_b_khigh_[a-z0-9_]+$','once')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
base=fullfile(root,'results','comparisons','two_dimensional_channel_closeout');
prior=fullfile(base,'20261023_s2c_a_1642');
formal=fullfile(base,'20261023_s2b_1305','formal_reference_v2');
out=fullfile(base,run_id); assert(~isfolder(out),'Refusing existing run');
freeze=jsondecode(fileread(fullfile(prior,'protocol','design_freeze.json')));
assert(strcmp(freeze.stage,'S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE'));
assert(abs(freeze.cost.charged_baseline_s-1408.660105)<1e-8);
assert(abs(freeze.cost.remaining_s-5791.339895)<1e-8);
formal_status=jsondecode(fileread(fullfile(formal,'summary.json')));
assert(strcmp(formal_status.status,'S2B_COMPLETE_REFERENCE_QUALIFIED'));
cfg_path=fullfile(prior,'protocol','K_HIGH_configuration.json');
cfg=jsondecode(fileread(cfg_path));
assert(strcmp(cfg.case_id,'K_HIGH') && cfg.A_m==.05 && cfg.K_rad_per_m==.2 && cfg.phi_rad==0);
hashes=readtable(fullfile(prior,'protocol','configuration_hashes.csv'),TextType='string');
assert(strcmp(closeout_s2b_hash(cfg_path),hashes.sha256(hashes.case_id=="K_HIGH")));
source_path=fullfile(prior,'protocol','frozen_source.mat');
line_path=fullfile(prior,'protocol','receivers_line.csv');
nine_path=fullfile(prior,'protocol','receivers_nine.csv');
assert(strcmp(closeout_s2b_hash(source_path),cfg.source_sha256));
assert(strcmp(closeout_s2b_hash(line_path),cfg.line_sha256));
assert(strcmp(closeout_s2b_hash(nine_path),cfg.nine_sha256));
source=load(source_path,'src');
formal_data=load(fullfile(formal,'reference_data.mat'),'U','x','z','ix9','src','flat_exact');
diagnostic=load(fullfile(prior,'mechanism_baseline','diagnostic_data.mat'),'plot_data');
assert(isequal(size(formal_data.U),[417 4 3]));
assert(max(abs(formal_data.U(:,3,:)-diagnostic.plot_data.U(:,3,:)),[],'all')==0);
assert(isequal(formal_data.src.coeff,source.src.coeff));
assert(all(formal_data.z(1:411)==3) && nnz(formal_data.z~=3)==6);
assert(numel(formal_data.ix9)==9);
line=readtable(line_path); nine=readtable(nine_path);
assert(height(line)==411 && height(nine)==9);
assert(max(abs(line.x_m-formal_data.x(1:411)))<1e-13);
assert(max(abs(nine.x_m-formal_data.x(formal_data.ix9)))<1e-13);
assert(max(abs(nine.z_m-formal_data.z(formal_data.ix9)))<1e-13);
gx=(-60:.01:60)';
[eta0,d10,d20]=local_geometry(gx,.05,.1,0);
[eta1,d11,d21]=local_geometry(gx,.05,.2,0);
saved=readtable(fullfile(prior,'protocol','candidate_geometry.csv'),TextType='string');
a=saved(saved.case_id=="BASE_M",:); b=saved(saved.case_id=="K_HIGH",:);
assert(height(a)==numel(gx) && height(b)==numel(gx));
assert(max(abs(a.x_m-gx))<1e-12 && max(abs(b.x_m-gx))<1e-12);
geometry_max=max(abs([eta0-a.eta_m;d10-a.slope;d20-a.curvature_per_m; ...
    eta1-b.eta_m;d11-b.slope;d21-b.curvature_per_m]));
assert(geometry_max<1e-12);
assert(all(eta1(abs(gx)>=50)==0) && all(d11(abs(gx)>=50)==0) && all(d21(abs(gx)>=50)==0));
% All five references are budgeted before any representative job.
nodes=[3840 3840 4480 3200 3736]';
levels=["L60_ppw12_qa96";"L60_ppw12_qa64";"L60_ppw14_qa96"; ...
    "L60_ppw10_qa96";"L70_ppw10_qa96"];
historical=[242.9793715;276.6899304;289.741466; ...
    242.9793715*(3200/3840)^2;242.9793715*(3736/3840)^2];
upper=1.2*historical;
estimated_increment=4*16*nodes.^2+128*1024^2;
bie_est=sum(upper); other_reserve=150; contingency=150;
total_est=bie_est+other_reserve+contingency;
cap=freeze.cost.suggested_S2c_B_cap_s;
projected_remaining=freeze.cost.remaining_s-total_est;
ram_pass=all(estimated_increment<=.7*available_physical_bytes);
cost_pass=total_est<=cap && projected_remaining>=freeze.cost.S3A_reserve_floor_s ...
    && all(upper<freeze.cost.single_job_cap_s);
mkdir(out); mkdir(fullfile(out,'protocol')); mkdir(fullfile(out,'P0')); mkdir(fullfile(out,'jobs'));
authorization=struct('goal','S2c-B K_HIGH only','authorized_case','K_HIGH', ...
    'A_m',.05,'K_rad_per_m',.2,'phi_rad',0,'explicit_user_authorization',true, ...
    'authorization_date','2026-10-23 Asia/Shanghai', ...
    'other_candidates_authorized',false,'S3_authorized',false, ...
    'frozen_configuration_sha256',closeout_s2b_hash(cfg_path), ...
    'upstream_design_freeze_sha256',closeout_s2b_hash(fullfile(prior,'protocol','design_freeze.json')));
closeout_s2b_json(fullfile(out,'protocol','authorization.json'),authorization);
if ram_pass && cost_pass, decision='P1_PASS_REPRESENTATIVE_ONLY';
else, decision='P1_BLOCKED_NO_PHYSICS'; end
gate=struct('decision',decision,'new_acoustic_seconds',0, ...
    'charged_prior_seconds',freeze.cost.charged_baseline_s,'core_remaining_s',freeze.cost.remaining_s, ...
    'subbudget_s',cap,'S3A_reserve_s',freeze.cost.S3A_reserve_floor_s, ...
    'BIE_five_upper_s',bie_est,'other_physics_reserve_s',other_reserve, ...
    'failure_contingency_s',contingency,'total_conservative_s',total_est, ...
    'projected_remaining_s',projected_remaining, ...
    'physical_available_bytes',available_physical_bytes, ...
    'measurement','Windows Memory Available Bytes immediately before P1', ...
    'max_increment_bytes',max(estimated_increment),'ram_reserve_fraction',.3, ...
    'RAM_pass',ram_pass,'cost_pass',cost_pass, ...
    'reference_levels',levels,'node_counts',nodes,'historical_stage_basis_s',historical, ...
    'upper_stage_s',upper,'single_job_cap_s',900, ...
    'next_step','Only L60 ppw12 qa96 representative; reestimate all remaining jobs from actual receipt');
closeout_s2b_json(fullfile(out,'protocol','P1_initial_gate.json'),gate);
column_names={'case_id','x_m','eta_m','slope','curvature_per_m'};
geometry=[table(repmat("BASE_M",numel(gx),1),gx,eta0,d10,d20,VariableNames=column_names); ...
    table(repmat("K_HIGH",numel(gx),1),gx,eta1,d11,d21,VariableNames=column_names)];
writetable(geometry,fullfile(out,'P0','geometry.csv'));
summary=table(["BASE_M";"K_HIGH"],[.1;.2],[2*pi/.1;2*pi/.2], ...
    [max(abs(eta0));max(abs(eta1))],[max(abs(d10));max(abs(d11))], ...
    [max(abs(d20));max(abs(d21))],[d10(gx==0);d11(gx==0)], ...
    'VariableNames',{'case_id','K_rad_per_m','period_m','max_abs_eta_m','max_abs_slope','max_abs_curvature_per_m','center_slope'});
writetable(summary,fullfile(out,'P0','geometry_summary.csv'));
Ubase=squeeze(formal_data.U(:,3,:));
Uflat=squeeze(formal_data.U(:,1,:));
raw=table((1:417)',formal_data.x,formal_data.z, ...
    real(Ubase(:,1)),imag(Ubase(:,1)),real(Ubase(:,2)),imag(Ubase(:,2)), ...
    real(Ubase(:,3)),imag(Ubase(:,3)), ...
    'VariableNames',{'sample_index','x_m','z_m','PE_real','PE_imag','BIE_real','BIE_imag','FullFK_real','FullFK_imag'});
writetable(raw,fullfile(out,'P0','BASE_M_fields.csv'));
save(fullfile(out,'P0','p0_data.mat'),'Ubase','Uflat','gx','eta0','d10','d20','eta1','d11','d21','-v7');
receipt=struct('status',decision,'new_physics',0, ...
    'BASE_M_source','S2b formal_reference_v2 and S2c-A diagnostic_data', ...
    'BASE_M_field_max_difference',0,'geometry_max_difference',geometry_max, ...
    'flat_outside_abs_x_50',true,'source_sha256',cfg.source_sha256, ...
    'S2b_reference_data_sha256',closeout_s2b_hash(fullfile(formal,'reference_data.mat')), ...
    'S2c_A_diagnostic_data_sha256',closeout_s2b_hash(fullfile(prior,'mechanism_baseline','diagnostic_data.mat')), ...
    'MATLAB_version',version,'configuration_sha256',authorization.frozen_configuration_sha256);
closeout_s2b_json(fullfile(out,'P0','receipt.json'),receipt);
p=struct('run_dir',out,'gate',gate,'receipt',receipt); disp(gate); disp(summary);
end

function [eta,d1,d2]=local_geometry(x,A,K,p)
r=abs(x); chi=double(r<=42); cp=zeros(size(x)); cpp=cp;
ii=r>42&r<50; t=(r(ii)-42)/8;
chi(ii)=1-10*t.^3+15*t.^4-6*t.^5;
cp(ii)=(-30*t.^2+60*t.^3-30*t.^4).*sign(x(ii))/8;
cpp(ii)=(-60*t+180*t.^2-120*t.^3)/64;
s=sin(K*x+p); c=cos(K*x+p);
eta=A*s.*chi; d1=A*(K*c.*chi+s.*cp);
d2=A*(-K^2*s.*chi+2*K*c.*cp+s.*cpp);
end
