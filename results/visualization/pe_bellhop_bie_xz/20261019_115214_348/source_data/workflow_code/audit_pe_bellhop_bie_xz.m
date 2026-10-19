function audit = audit_pe_bellhop_bie_xz(run_dir, rebuilt_figure_dir)
%AUDIT_PE_BELLHOP_BIE_XZ Check raw-field CSV, transforms, source and rebuild.
t=load(fullfile(run_dir,'plot_data','xz_data.mat'),'data'); d=t.data;
assert(numel(d.x_m)==411 && numel(d.z_m)==45 && all(diff(d.x_m)>0) && all(diff(d.z_m)>0));
manifest=readtable(fullfile(run_dir,'source_data','input_manifest.csv'),'TextType','string');
for j=1:height(manifest)
    assert(strcmpi(file_hash(fullfile(run_dir,'source_data',manifest.snapshot(j))),manifest.sha256(j)));
end
T=readtable(fullfile(run_dir,'plot_data','xz_fields.csv')); [X,Z]=ndgrid(d.x_m,d.z_m);
assert(max(abs(T.x_m-X(:)))<1e-12 && max(abs(T.z_m-Z(:)))<1e-12);
max_csv=0;
raw_pe=complex(T.PE_envelope_real,T.PE_envelope_imag);
raw_bh=complex(T.BH_raw_SHD_real,T.BH_raw_SHD_imag);
assert(max(abs(raw_pe-d.PE_envelope(:)))<1e-12 && max(abs(raw_bh-d.BH_raw_SHD(:)))<1e-12);
for j=1:3
    method=d.methods{j}; u=d.(method); value=complex(T.([method '_real']),T.([method '_imag']));
    assert(isequal(size(u),[411 45]) && all(isfinite(u),'all'));
    max_csv=max(max_csv,max(abs(u(:)-value)));
    assert(isequal(abs(u),d.magnitude.(method)) && isequal(angle(u),d.phase_rad.(method)));
end
assert(max_csv<1e-12);
k=2*pi*4000/1500;
carrier_error=norm(d.PE_carrier_restored-d.PE_envelope.*exp(1i*k*(100+d.z_m.')),'fro');
assert(carrier_error==0 && norm(d.BH_native-conj(d.BH_raw_SHD),'fro')==0);
assert(~d.source_diagnostic.fitted_scalar_applied && ~d.source_diagnostic.phase_alignment_applied);
trans=load(fullfile(run_dir,'figures','plot_transforms.mat'),'transforms'); tr=trans.transforms;
assert(~tr.fields_use_ratio && isequal(tr.phase_limits_rad,[-pi pi]));
assert(tr.ceilings(1)==tr.ceilings(3)); max_rgb_error=0;
for j=1:3
    method=d.methods{j}; u=d.(method);
    expected=hsv2rgb(cat(3,(angle(u)+pi)/(2*pi),ones(size(u)),min(abs(u)/tr.ceilings(j),1)));
    max_rgb_error=max(max_rgb_error,max(abs(expected-tr.complex_rgb.(method)),[],'all'));
end
assert(max_rgb_error==0);
names={'XZ01_raw_magnitude','XZ02_raw_phase','XZ03_complex_field','XZ04_raw_profiles'};
same=false(1,4);
for j=1:4
    a=imread(fullfile(run_dir,'figures',[names{j} '.png']));
    b=imread(fullfile(rebuilt_figure_dir,[names{j} '.png']));
    assert(size(a,1)==1350 && size(a,2)==2400); same(j)=isequal(a,b);
end
assert(all(same));
assert(d.new_Bellhop_jobs==1 && d.new_BIE_boundary_solves==0);
audit=struct('passed',true,'CSV_max_abs_difference',max_csv,'carrier_formula_difference',carrier_error, ...
    'RGB_formula_difference',max_rgb_error,'identical_PNG_rebuilds',nnz(same), ...
    'grid_points_per_method',numel(X),'new_Bellhop_jobs',1,'new_BIE_boundary_solves',0, ...
    'BIE_z3_closure',d.checks.BIE_z3_closure,'PE_z3_closure',d.checks.PE_z3_closure, ...
    'BH_z3_closure',d.checks.BH_z3_closure,'matlab_version',version);
fid=fopen(fullfile(run_dir,'metrics','artifact_audit.json'),'w'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); disp(audit);
end

function value=file_hash(path)
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid)); md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), md.update(fread(fid,1048576,'*uint8')); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
