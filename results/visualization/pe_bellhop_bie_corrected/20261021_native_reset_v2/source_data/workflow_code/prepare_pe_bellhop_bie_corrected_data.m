function output = prepare_pe_bellhop_bie_corrected_data(run_id)
%PREPARE_PE_BELLHOP_BIE_CORRECTED_DATA Reset fixed-case visualization data.
% Validation-only saved-start propagation; no source, BIE or Bellhop solve.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1, run_id=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')); end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_-]+$','once')),'Invalid run ID.');
out=fullfile(root,'results','visualization','pe_bellhop_bie_corrected',run_id);
assert(~isfolder(out),'Refusing to overwrite a completed or partial run.');
base=fullfile(root,'results','visualization');
line_root=fullfile(base,'pe_bellhop_bie_comparison','20261018_194331_523');
xz_root=fullfile(base,'pe_bellhop_bie_xz','20261019_115214_348');
sources={fullfile(line_root,'plot_data','comparison_data.mat'), ...
    fullfile(xz_root,'plot_data','xz_data.mat'), ...
    fullfile(xz_root,'source_data','numeric_inputs.mat')};
assert(all(cellfun(@isfile,sources)),'Required saved numerical source is missing.');
s=load(sources{1},'data'); old=s.data;
s=load(sources{2},'data'); old_xz=s.data;
s=load(sources{3},'inputs'); inputs=s.inputs;
assert(isequal(old.x_m(:),inputs.xfull_m(:)));
assert(isequal(old_xz.x_m(:),old.x_m(old.mask)));
assert(isequal(old_xz.z_m(:),(.5:.125:6).'));
assert(isequal(logical(inputs.original_M99(:)),logical(old.mask(:))));
k=2*pi*4000/1500; kx=inputs.pe.kx_rad_per_m(:);
kz=sqrt(complex(k^2-kx.^2,0)); inc=inputs.pe.surface_incident_field(:);
assert(numel(inc)==numel(old.x_m) && numel(kx)==numel(inc));
H3=exp(1i*3*(kz-k)); flat=ifft(fft(-inc).*H3);
flat_closure=norm(flat-inputs.PE_flat_envelope(:))/norm(flat);
assert(flat_closure<1e-11,'Flat PE control has changed.');
data=old; data.schema_version='corrected-native-v1';
data.convention=['All G in native exp(-i omega t) representation; physical ' ...
    'surface coordinate z=eta positive downward; corrected PE screen exp(-2iketa).'];
data.legacy_metrics=old.metrics;
data.metrics=old.metrics;
data.metrics.complex_l2(:)=NaN;
data.metrics.phase_rms_rad(:)=NaN;
data.metrics.magnitude_ratio_rms_db(:)=NaN;
data.cases=old.cases;
[data.cases.G_PE_legacy]=deal([]);
[data.cases.G_BIE_legacy]=deal([]);
[data.cases.G_BH_legacy]=deal([]);
data.correction=struct('origin','saved incident field and fixed old geometry', ...
    'old_PE_screen','exp(+2iketa)','new_PE_screen','exp(-2iketa)', ...
    'physical_surface','z=eta, depth positive downward', ...
    'time_convention','exp(-i omega t)', ...
    'flat_relative_l2',flat_closure,'new_source_propagations',0, ...
    'new_BIE_solves',0,'new_BIE_field_evaluations',0,'new_Bellhop_jobs',0);
comparison_rows=cell(9,1); decomposition_rows=cell(9,1); row=0;
for j=1:3
    c=data.cases(j); eta=c.eta_m(:); mask=c.plot_mask(:);
    assert(isequal(mask,old.mask(:)) && all(abs(old.x_m(mask))<42));
    plus=ifft(fft(-inc.*exp(1i*2*k*eta)).*H3)./flat;
    assert(max(abs(plus(mask)-old.cases(j).G_PE(mask)))<1e-11, ...
        'Old PE result does not reproduce.');
    rough_spec=fft(-inc.*exp(-1i*2*k*eta));
    minus=ifft(rough_spec.*H3)./flat;
    c.G_PE_legacy=old.cases(j).G_PE;
    c.G_BIE_legacy=old.cases(j).G_BIE;
    c.G_BH_legacy=old.cases(j).G_BH;
    c.G_PE=minus; c.G_BIE=c.G_BIE_native; c.G_BH=c.G_BH_native;
    c.G_BH_5001=conj(c.G_BH_5001);
    assert(all(isfinite(real(c.G_PE(mask))) & isfinite(imag(c.G_PE(mask)))));
    assert(all(isfinite(real(c.G_BIE(mask))) & isfinite(imag(c.G_BIE(mask)))));
    assert(all(isfinite(real(c.G_BH(mask))) & isfinite(imag(c.G_BH(mask)))));
    pair_fields={'G_PE','G_BH','G_PE'}; ref_fields={'G_BIE','G_BIE','G_BH'};
    pair_names={'PE_BIE','BH_BIE','PE_BH'}; w=data.weights(mask); w=w/sum(w);
    for p=1:3
        row=row+1; a=c.(pair_fields{p}); b=c.(ref_fields{p});
        residual=struct('pair',pair_names{p},'complex_abs',abs(a-b), ...
            'phase_rad',angle(a.*conj(b)), ...
            'magnitude_ratio_db',20*log10(abs(a)./abs(b)));
        c.residuals(p)=residual;
        aa=a(mask); bb=b(mask); dr=residual.phase_rad(mask);
        dm=abs(aa)-abs(bb);
        complex_l2=sqrt(sum(w.*abs(aa-bb).^2)/sum(w.*abs(bb).^2));
        mag_rel=sqrt(sum(w.*dm.^2)/sum(w.*abs(bb).^2));
        phase_rms=sqrt(sum(w.*dr.^2));
        db_rms=sqrt(sum(w.*residual.magnitude_ratio_db(mask).^2));
        data.metrics.complex_l2(row)=complex_l2;
        data.metrics.phase_rms_rad(row)=phase_rms;
        data.metrics.magnitude_ratio_rms_db(row)=db_rms;
        data.metrics.sample_count(row)=nnz(mask);
        comparison_rows{row}=table(string(c.name),string(pair_names{p}), ...
            complex_l2,mag_rel,phase_rms,db_rms, ...
            old.metrics.complex_l2(row),old.metrics.phase_rms_rad(row), ...
            old.metrics.magnitude_ratio_rms_db(row), ...
            'VariableNames',{'case_name','pair','complex_l2','magnitude_relative_l2', ...
            'phase_rms_rad','magnitude_ratio_rms_db','legacy_complex_l2', ...
            'legacy_phase_rms_rad','legacy_magnitude_ratio_rms_db'});
        % Exact pointwise split of |A-B|^2 into magnitude and phase terms.
        den=sum(w.*abs(bb).^2);
        mag_term=sum(w.*dm.^2)/den;
        phase_term=sum(w.*(2*abs(aa).*abs(bb).*(1-cos(dr))))/den;
        assert(abs(complex_l2^2-mag_term-phase_term)<1e-12);
        decomposition_rows{row}=table(string(c.name),string(pair_names{p}), ...
            mag_term,phase_term,mag_term/(mag_term+phase_term), ...
            phase_term/(mag_term+phase_term), ...
            'VariableNames',{'case_name','pair','magnitude_term','phase_term', ...
            'magnitude_fraction','phase_fraction'});
    end
    c.transforms.PE=transform_curve(c.G_PE);
    c.transforms.BIE=transform_curve(c.G_BIE);
    c.transforms.BH=transform_curve(c.G_BH);
    data.cases(j)=c;
end
comparison=vertcat(comparison_rows{:}); decomposition=vertcat(decomposition_rows{:});
assert(max(abs(data.metrics.complex_l2-comparison.complex_l2))<1e-14);
xz=old_xz; strong=data.cases(3); n=numel(old_xz.z_m);
PE_envelope=complex(zeros(numel(old_xz.x_m),n));
for iz=1:n
    H=exp(1i*old_xz.z_m(iz)*(kz-k)); u=ifft(rough_spec.*H);
    PE_envelope(:,iz)=u(old.mask);
end
xz.PE_envelope=PE_envelope;
xz.PE_carrier_restored=PE_envelope.*exp(1i*k*(100+old_xz.z_m(:).'));
xz.magnitude.PE_carrier_restored=abs(xz.PE_carrier_restored);
xz.phase_rad.PE_carrier_restored=angle(xz.PE_carrier_restored);
xz.checks.PE_z3_closure=norm(PE_envelope(:,old_xz.z_m==3)- ...
    strong.G_PE(old.mask).*flat(old.mask))/norm(PE_envelope(:,old_xz.z_m==3));
assert(xz.checks.PE_z3_closure<1e-11);
xz.conventions.PE=['Corrected native field: exp(-2iketa) screen, FFT to z, ' ...
    'carrier restored exp(+ik(100+z))'];
xz.schema_version='corrected-native-v1'; xz.new_PE_full_runs=0;
xz.PE_saved_start_depth_evaluations=n;
assert(isequal(size(xz.PE_envelope),[411 45]));
output=struct('run_dir',out,'line_source',sources{1},'xz_source',sources{2}, ...
    'incident_source',sources{3},'flat_relative_l2',flat_closure, ...
    'PE_z3_closure',xz.checks.PE_z3_closure,'new_saved_start_FFT_fields',3+n, ...
    'new_source_propagations',0,'new_BIE_solves',0,'new_Bellhop_jobs',0);
mkdir(out); mkdir(fullfile(out,'plot_data')); mkdir(fullfile(out,'source_data'));
mkdir(fullfile(out,'metrics'));
workflow=fullfile(out,'source_data','workflow_code'); mkdir(workflow);
input_dir=fullfile(out,'source_data','input_mat'); mkdir(input_dir);
packaged={'old_comparison_data.mat','old_xz_data.mat','numeric_inputs.mat'};
for i=1:3, copyfile(sources{i},fullfile(input_dir,packaged{i})); end
renderers={'generate_pe_bellhop_bie_corrected_visuals.m', ...
    'generate_pe_bellhop_bie_comparison_visuals.m', ...
    'generate_pe_bellhop_bie_xz_visuals.m'};
for i=1:numel(renderers)
    source=fullfile(root,'scripts','reporting',renderers{i});
    assert(isfile(source)); copyfile(source,fullfile(workflow,renderers{i}));
end
validation_code={'prepare_pe_bellhop_bie_corrected_data.m', ...
    'audit_pe_bellhop_bie_corrected.m'};
for i=1:numel(validation_code)
    source=fullfile(root,'scripts','validation',validation_code{i});
    assert(isfile(source)); copyfile(source,fullfile(workflow,validation_code{i}));
end
save(fullfile(out,'plot_data','comparison_data.mat'),'data','-v7.3');
for j=1:3
    c=data.cases(j); x_m=data.x_m(:); eta_m=c.eta_m(:);
    fixed_M99=logical(c.plot_mask(:));
    G_PE_real=real(c.G_PE(:)); G_PE_imag=imag(c.G_PE(:));
    G_BIE_real=real(c.G_BIE(:)); G_BIE_imag=imag(c.G_BIE(:));
    G_BH_real=real(c.G_BH(:)); G_BH_imag=imag(c.G_BH(:));
    legacy_PE_real=real(c.G_PE_legacy(:));
    legacy_PE_imag=imag(c.G_PE_legacy(:));
    T=table(x_m,eta_m,fixed_M99,G_PE_real,G_PE_imag, ...
        G_BIE_real,G_BIE_imag,G_BH_real,G_BH_imag, ...
        legacy_PE_real,legacy_PE_imag);
    writetable(T,fullfile(out,'plot_data',[c.name '_line.csv']));
end
data=xz; save(fullfile(out,'plot_data','xz_data.mat'),'data','-v7.3');
x_m=repmat(xz.x_m(:),numel(xz.z_m),1);
z_m=kron(xz.z_m(:),ones(numel(xz.x_m),1));
PE_real=real(xz.PE_carrier_restored(:)); PE_imag=imag(xz.PE_carrier_restored(:));
BIE_real=real(xz.BIE_native(:)); BIE_imag=imag(xz.BIE_native(:));
BH_real=real(xz.BH_native(:)); BH_imag=imag(xz.BH_native(:));
writetable(table(x_m,z_m,PE_real,PE_imag,BIE_real,BIE_imag, ...
    BH_real,BH_imag),fullfile(out,'plot_data','xz_fields.csv'));
writetable(comparison,fullfile(out,'metrics','corrected_vs_legacy.csv'));
writetable(decomposition,fullfile(out,'metrics','complex_error_decomposition.csv'));
% Store immutable numerical inputs by hash/path; redraw only reads plot_data.
relative_paths=string(fullfile('source_data','input_mat',packaged(:)));
source_table=table(string(sources(:)),relative_paths,strings(3,1), ...
    'VariableNames',{'source_path','packaged_path','SHA256'});
for i=1:3, source_table.SHA256(i)=string(file_sha256(sources{i})); end
writetable(source_table,fullfile(out,'source_data','source_hashes.csv'));
fid=fopen(fullfile(out,'prepare_audit.json'),'w','n','UTF-8');
assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); clear cl
disp(comparison); disp(output);
end

function t=transform_curve(g)
t=struct('magnitude',abs(g),'phase_rad',angle(g), ...
    'phase_mrad',angle(g)*1e3);
gap=false(size(g)); ph=t.phase_rad; gap(2:end)=abs(diff(ph))>pi;
shown=t.phase_mrad; shown(gap)=NaN;
t.phase_wrap_gap=gap; t.plotted_phase_mrad=shown;
end

function h=file_sha256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8'); md.update(bytes);
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
