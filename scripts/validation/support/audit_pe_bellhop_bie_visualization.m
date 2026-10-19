function audit = audit_pe_bellhop_bie_visualization(run_dir, rebuild_dir)
%AUDIT_PE_BELLHOP_BIE_VISUALIZATION Independent saved-artifact/CSV/pixel checks.
q=load(fullfile(run_dir,'plot_data','comparison_data.mat'),'data'); d=q.data;
checks=struct('all_jobs',true,'source_snapshots',true,'native_mapping',true, ...
    'shared_mask',true,'csv_roundtrip',true,'rebuild_pixels',true,'binary_hash',true);
manifest=readtable(fullfile(run_dir,'source_data','input_manifest.csv'),'TextType','string');
for j=1:height(manifest)
    path=fullfile(run_dir,'source_data',manifest.snapshot(j));
    assert(strcmp(sha256(path),manifest.sha256(j)),'Source snapshot hash mismatch.');
end
origin=jsondecode(fileread(fullfile(run_dir,'source_data','binary_origin.json')));
assert(strcmp(sha256(fullfile(run_dir,'source_data','bellhop_internal_wall.exe')),origin.binary_sha256));
max_csv_difference=0; max_mapping_difference=0;
for ci=1:3
    c=d.cases(ci); assert(isequal(c.plot_mask,d.mask));
    T=readtable(fullfile(run_dir,'plot_data',[c.name '.csv']));
    observed=[complex(T.PE_real,T.PE_imag),complex(T.BH_real,T.BH_imag),complex(T.BIE_real,T.BIE_imag)];
    expected=[c.G_PE,c.G_BH,c.G_BIE];
    good=c.valid_mask; assert(isequal(isnan(observed),isnan(expected)));
    max_csv_difference=max(max_csv_difference,max(abs(observed(good,:)-expected(good,:)),[],'all'));
    assert(max(abs(T.x_m-d.x_m))<1e-12 && isequal(logical(T.M99_mask),c.plot_mask));
    assert(max(abs(T.incident_weight-d.weights))<1e-12);
    assert(all(isfinite(observed(good,:)),'all') && all(good(d.mask)));
    assert(isequal(logical(T.ratio_valid),good));
    for method={'PE','BH','BIE'}
        tag=method{1}; z=c.(['G_' tag]); tr=c.transforms.(tag);
        assert(isequal(tr.phase_wrap_gap,[false;abs(diff(angle(z)))>pi]));
        finite=isfinite(z);
        assert(max(abs(tr.magnitude(finite)-abs(z(finite))))==0 ...
            && max(abs(tr.phase_rad(finite)-angle(z(finite))))==0);
    end
    low=load(fullfile(run_dir,'source_data',sprintf('job_%02d.mat',ci*2+1)),'job_result');
    hi=load(fullfile(run_dir,'source_data',sprintf('job_%02d.mat',ci*2+2)),'job_result');
    fl=load(fullfile(run_dir,'source_data','job_01.mat'),'job_result');
    fh=load(fullfile(run_dir,'source_data','job_02.mat'),'job_result');
    mapped=conj(hi.job_result.native_field(good)./fh.job_result.native_field(good));
    low_good=c.BH5001_valid_mask;
    mapped_low=conj(low.job_result.native_field(low_good)./fl.job_result.native_field(low_good));
    max_mapping_difference=max([max_mapping_difference;abs(mapped-c.G_BH(good));abs(mapped_low-c.G_BH_5001(low_good))]);
end
assert(max_csv_difference<1e-12 && max_mapping_difference<1e-14);
for j=1:8
    q=load(fullfile(run_dir,'source_data',sprintf('job_%02d.mat',j)),'job_result');
    assert(q.job_result.geometry.pass && all(isfinite(q.job_result.raw_SHD)));
    assert(norm(q.job_result.native_field-conj(q.job_result.raw_SHD))==0);
end
list={'01_environment','bellhop_bie/BB01_fields','bellhop_bie/BB02_residuals', ...
    'bellhop_bie/BB03_metrics','three_way/TR01_fields','three_way/TR02_residuals','three_way/TR03_metrics'};
pixels_equal=false(7,1);
for j=1:7
    a=imread(fullfile(run_dir,'figures',[list{j} '.png']));
    b=imread(fullfile(rebuild_dir,[list{j} '.png']));
    assert(size(a,1)==1350 && size(a,2)==2400); pixels_equal(j)=isequal(a,b);
end
assert(all(pixels_equal),'Rendered PNG pixels differ on rebuild.');
audit=struct('checks',checks,'passed',all(structfun(@(x)x,checks)), ...
    'csv_max_abs_difference',max_csv_difference,'native_mapping_max_abs_difference',max_mapping_difference, ...
    'identical_PNG_rebuilds',nnz(pixels_equal),'new_PE_solves',0,'new_BIE_solves',0, ...
    'new_Bellhop_jobs',8,'matlab_version',version);
fid=fopen(fullfile(run_dir,'metrics','artifact_audit.json'),'w'); clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); disp(audit);
end

function value=sha256(path)
fid=fopen(path,'rb'); assert(fid>=0); clean=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), md.update(fread(fid,1048576,'*uint8')); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
