function audit=audit_pe_bellhop_bie_corrected(run_dir,rebuild_dir)
%AUDIT_PE_BELLHOP_BIE_CORRECTED Independent saved-data and redraw check.
if nargin<2, rebuild_dir=fullfile(run_dir,'rebuilt_figures'); end
s=load(fullfile(run_dir,'plot_data','comparison_data.mat'),'data'); d=s.data;
s=load(fullfile(run_dir,'plot_data','xz_data.mat'),'data'); xz=s.data;
assert(strcmp(d.schema_version,'corrected-native-v1'));
assert(strcmp(xz.schema_version,'corrected-native-v1'));
M=readtable(fullfile(run_dir,'metrics','corrected_vs_legacy.csv'),'Delimiter',',');
D=readtable(fullfile(run_dir,'metrics','complex_error_decomposition.csv'),'Delimiter',',');
assert(height(M)==9 && height(D)==9);
line_max=0; metric_max=0; split_max=0; n=0;
for j=1:3
    c=d.cases(j); mask=c.plot_mask(:);
    assert(isequal(mask(:),logical(d.mask(:))) && nnz(mask)==411);
    T=readtable(fullfile(run_dir,'plot_data',[c.name '_line.csv']), ...
        'Delimiter',',');
    assert(max(abs(T.x_m-d.x_m(:)))<1e-12);
    assert(isequal(logical(T.fixed_M99),mask));
    csv_fields={complex(T.G_PE_real,T.G_PE_imag), ...
        complex(T.G_BIE_real,T.G_BIE_imag), ...
        complex(T.G_BH_real,T.G_BH_imag), ...
        complex(T.legacy_PE_real,T.legacy_PE_imag)};
    mat_fields={c.G_PE,c.G_BIE,c.G_BH,c.G_PE_legacy};
    for p=1:4
        line_max=max(line_max,max(abs(csv_fields{p}(mask)-mat_fields{p}(mask))));
    end
    A={c.G_PE,c.G_BH,c.G_PE}; B={c.G_BIE,c.G_BIE,c.G_BH};
    w=d.weights(mask); w=w/sum(w);
    for p=1:3
        n=n+1; a=A{p}(mask); b=B{p}(mask); phi=angle(a.*conj(b));
        den=sum(w.*abs(b).^2);
        E=sqrt(sum(w.*abs(a-b).^2)/den);
        E_mag=sqrt(sum(w.*(abs(a)-abs(b)).^2)/den);
        E_phase=sqrt(sum(w.*phi.^2));
        E_db=sqrt(sum(w.*(20*log10(abs(a)./abs(b))).^2));
        metric_max=max(metric_max,max(abs([E-M.complex_l2(n), ...
            E_mag-M.magnitude_relative_l2(n), ...
            E_phase-M.phase_rms_rad(n),E_db-M.magnitude_ratio_rms_db(n)])));
        tmag=sum(w.*(abs(a)-abs(b)).^2)/den;
        tph=sum(w.*2.*abs(a).*abs(b).*(1-cos(phi)))/den;
        split_max=max(split_max,max(abs([tmag-D.magnitude_term(n), ...
            tph-D.phase_term(n),E^2-tmag-tph])));
    end
end
assert(line_max<1e-12 && metric_max<1e-12 && split_max<1e-12);
T=readtable(fullfile(run_dir,'plot_data','xz_fields.csv'),'Delimiter',',');
assert(height(T)==411*45);
assert(max(abs(T.x_m-repmat(xz.x_m(:),45,1)))<1e-12);
assert(max(abs(T.z_m-kron(xz.z_m(:),ones(411,1))))<1e-12);
xz_fields={complex(T.PE_real,T.PE_imag), ...
    complex(T.BIE_real,T.BIE_imag),complex(T.BH_real,T.BH_imag)};
mat_fields={xz.PE_carrier_restored(:),xz.BIE_native(:),xz.BH_native(:)};
xz_max=0;
for p=1:3
    xz_max=max(xz_max,max(abs(xz_fields{p}-mat_fields{p})));
end
assert(xz_max<1e-12);
zi=find(xz.z_m==3); assert(isscalar(zi));
carrier=exp(1i*(2*pi*4000/1500)*(100+3));
denom=norm(xz.PE_carrier_restored(:,zi));
slice_closure=norm(xz.PE_carrier_restored(:,zi)- ...
    xz.PE_envelope(:,zi)*carrier)/denom;
assert(slice_closure<1e-12);
files=readtable(fullfile(run_dir,'figures','image_metadata.csv'), ...
    'Delimiter',',');
assert(height(files)==12);
pixel_difference=0;
for j=1:height(files)
    name=char(string(files.name(j)));
    a=imread(fullfile(run_dir,'figures',[name '.png']));
    b=imread(fullfile(rebuild_dir,[name '.png']));
    assert(isequal(size(a),[1350 2400 3]));
    pixel_difference=pixel_difference+nnz(a~=b);
    assert(isfile(fullfile(run_dir,'figures',[name '.pdf'])));
    assert(isfile(fullfile(rebuild_dir,[name '.pdf'])));
end
assert(pixel_difference==0);
H=readtable(fullfile(run_dir,'source_data','source_hashes.csv'), ...
    'Delimiter',',');
assert(height(H)==3);
source_hash_matches=0;
for j=1:height(H)
    p=fullfile(run_dir,char(string(H.packaged_path(j))));
    source_hash_matches=source_hash_matches+strcmpi( ...
        file_sha256(p),char(string(H.SHA256(j))));
end
assert(source_hash_matches==3);
audit=struct('status','PASS_SAVED_DATA_AND_PIXEL_REDRAW', ...
    'run_dir',run_dir,'line_csv_max_abs',line_max,'metric_max_abs',metric_max, ...
    'decomposition_max_abs',split_max,'xz_csv_max_abs',xz_max, ...
    'PE_carrier_slice_closure',slice_closure, ...
    'redraw_pixel_differences',pixel_difference, ...
    'source_hash_matches',source_hash_matches, ...
    'PNG_PDF_pairs',height(files),'new_solver_runs',0);
fid=fopen(fullfile(run_dir,'metrics','corrected_audit.json'),'w','n','UTF-8');
assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); clear cl
disp(audit);
end

function h=file_sha256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8'); md.update(bytes);
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
