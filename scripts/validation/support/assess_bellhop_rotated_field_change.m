function receipt = assess_bellhop_rotated_field_change(coarse_file,fine_file,output_dir,options)
%ASSESS_BELLHOP_ROTATED_FIELD_CHANGE Stratified complex/phase convergence.
% Normalize point errors by sum of component magnitudes, including at nulls;
% phase gate uses |p|>=.05*envelope. No fitted phase, gain, or interpolation.
arguments
    coarse_file (1,:) char
    fine_file (1,:) char
    output_dir (1,:) char
    options.label (1,:) char='numerical_change'
    options.relative_tolerance (1,1) double {mustBePositive}=.02
    options.envelope_tolerance (1,1) double {mustBePositive}=.02
    options.phase_tolerance_rad (1,1) double {mustBePositive}=.05
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
a=load(coarse_file,'total','plan'); b=load(fine_file,'total','envelope','plan');
assert(isequal(a.plan.depths_m,b.plan.depths_m) && isequal(a.plan.x_m,b.plan.x_m),'Receiver grids differ.');
valid=~b.plan.source_singularity_mask;
assert(all(isfinite(a.total(valid))) && all(isfinite(b.total(valid))) && all(b.envelope(valid)>0));
nz=numel(b.plan.depths_m); nx=numel(b.plan.x_m);
z=b.plan.depths_m(:)+zeros(1,nx); x=zeros(nz,1)+b.plan.x_m;
r=hypot(x,z-b.plan.physical_contract.source_z_m);
regions={valid,valid & z<=4900,valid & z>=4900 & r>=1,valid & x==0 & r>=1,valid & r<1};
names={'all','bulk_z0_4900','near_bottom_r_ge1','axis_r_ge1','source_vicinity_r_lt1'};
metrics=table(strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),false(0,1), ...
    'VariableNames',{'region','points','relative_L2','max_envelope_normalized','max_phase_rad','pass'});
for j=1:numel(regions)
    use=regions{j}; if ~any(use,'all'), continue; end
    delta=a.total(use)-b.total(use);
    relative=norm(delta)/norm(b.total(use)); maximum=max(abs(delta)./b.envelope(use));
    phase_use=use & abs(b.total)>=.05*b.envelope;
    ph=max(abs(angle(a.total(phase_use).*conj(b.total(phase_use)))));
    if isempty(ph), ph=0; end
    pass=relative<=options.relative_tolerance && maximum<=options.envelope_tolerance && ph<=options.phase_tolerance_rad;
    metrics(end+1,:)={names{j},nnz(use),relative,maximum,ph,pass}; %#ok<AGROW>
end
receipt=struct('status','COMPLETE','label',options.label,'failed_checks',nnz(~metrics.pass), ...
    'elapsed_s',toc(t),'result_dir',output_dir,'settings',options, ...
    'coarse_file',coarse_file,'coarse_sha256',black_file_hash(coarse_file), ...
    'fine_file',fine_file,'fine_sha256',black_file_hash(fine_file));
writetable(metrics,fullfile(output_dir,'checks.csv')); save(fullfile(output_dir,'result.mat'),'receipt','metrics');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('CONVERGENCE %s failed_regions=%d max_envelope=%.9g\n',options.label,receipt.failed_checks,max(metrics.max_envelope_normalized));
end