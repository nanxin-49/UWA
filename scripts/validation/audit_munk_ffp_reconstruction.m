function receipt=audit_munk_ffp_reconstruction(run_id)
%AUDIT_MUNK_FFP_RECONSTRUCTION Replay named raw samples using frozen source.
arguments
    run_id (1,:) char = '20261112_full_01'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(fullfile(root,'scripts','validation','support'),'-begin');
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id);
assert(~isfile(fullfile(d,'reconstruction_receipt.json'))); [~,mem]=memory;
assert(mem.PhysicalMemory.Available>=2*2^30,'Reconstruction block requires 2GiB start RAM');
t=tic; diary(fullfile(d,'reconstruction_diary.log'));
receipt=struct('stage','saved_field_reconstruction','status','RUNNING','elapsed_s',0,'exit_status',NaN);
try
    s=load(fullfile(d,'J07.mat'),'output'); o=s.output; c=o.cfg;
    depths=[0 700 1900 2950 4500 4949.97 4950 4950.03 4999.99 5000]; ranges=[-5 -1.15 0 .35 5];
    [distz,iz]=min(abs(o.depth_m-depths),[],1); [distr,ix]=min(abs(o.range_m-ranges.'),[],2);
    assert(max(distz)<1e-8 && max(distr)<1e-8); c.depth_m=o.depth_m(iz); c.range_m=o.range_m(ix);
    truth=o.green(iz,ix); clear s o;
    reconstructed=solve_layered_fluid_ffp(c); valid=isfinite(truth);
    rel=norm(reconstructed.green(valid)-truth(valid))/norm(truth(valid)); absmax=max(abs(reconstructed.green(valid)-truth(valid)));
    per_depth=zeros(numel(iz),1);
    for j=1:numel(iz)
        good=isfinite(truth(j,:)); den=norm(truth(j,good));
        if den>0, per_depth(j)=norm(reconstructed.green(j,good)-truth(j,good))/den;
        else, per_depth(j)=max(abs(reconstructed.green(j,good))); end
    end
    assert(max(per_depth)<=1e-10,'Current-source replay differs from raw saved field');
    save(fullfile(d,'reconstructed_samples.mat'),'reconstructed','truth','per_depth','iz','ix','-v7');
    receipt.status='CURRENT_SOURCE_RECONSTRUCTION_PASS'; receipt.exit_status=0;
    receipt.metrics=struct('sample_count',numel(truth),'relative_l2',rel,'max_absolute',absmax,'max_per_depth_relative',max(per_depth));
catch err
    receipt.status='RECONSTRUCTION_FAILED'; receipt.exit_status=1; receipt.error=getReport(err,'extended','hyperlinks','off');
end
receipt.elapsed_s=toc(t); fid=fopen(fullfile(d,'reconstruction_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,'PrettyPrint',true)); fclose(fid); diary off;
fprintf('%s\n',jsonencode(receipt));
end
