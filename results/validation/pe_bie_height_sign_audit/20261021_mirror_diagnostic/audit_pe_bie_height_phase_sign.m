function audit = audit_pe_bie_height_phase_sign(run_id)
%AUDIT_PE_BIE_HEIGHT_PHASE_SIGN Bounded saved-start sign diagnostic.
% Reevaluate the same 3 m FFT operator with +/- height screens. No source
% propagation, BIE solve/evaluation, Bellhop job or production edit is made.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1, run_id=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')); end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_-]+$','once')),'Invalid run ID.');
out=fullfile(root,'results','validation','pe_bie_height_sign_audit',run_id);
assert(~isfolder(out),'Refusing to overwrite an existing audit.');
base=fullfile(root,'results','visualization');
line_file=fullfile(base,'pe_bellhop_bie_comparison','20261018_194331_523', ...
    'plot_data','comparison_data.mat');
input_file=fullfile(base,'pe_bellhop_bie_xz','20261019_115214_348', ...
    'source_data','numeric_inputs.mat');
t=load(line_file,'data'); d=t.data;
t=load(input_file,'inputs'); in=t.inputs;
k=2*pi*4000/1500; z=3;
kz=sqrt(complex(k^2-in.pe.kx_rad_per_m(:).^2,0));
H=exp(1i*z*(kz-k)); inc=in.pe.surface_incident_field(:);
flat=ifft(fft(-inc).*H);
flat_closure=norm(flat-in.PE_flat_envelope)/norm(flat);
assert(flat_closure<1e-11);
assert(isequal(d.x_m(:),in.xfull_m(:)));
amplitude=abs(inc);
amplitude_prime=real(ifft(1i*in.pe.kx_rad_per_m(:).*fft(amplitude)));
cases=struct([]);
row_template=metrics([1;2],[1;2],[.5;.5]);
row_template.case_name=''; row_template.comparison='';
rows=repmat(row_template,9,1);
for j=1:3
    c=d.cases(j); eta=c.eta_m(:); mask=c.plot_mask;
    w=d.weights(mask); w=w/sum(w);
    plus=ifft(fft(-inc.*exp(1i*2*k*eta)).*H)./flat;
    minus=ifft(fft(-inc.*exp(-1i*2*k*eta)).*H)./flat;
    saved_closure=max(abs(plus(mask)-c.G_PE(mask)));
    assert(saved_closure<1e-11);
    native=c.G_BIE_native;
    assert(max(abs(native(mask)-conj(c.G_BIE(mask))))<1e-12);
    fields={plus,plus,minus}; refs={conj(native),native,native};
    labels={'legacy_plus_vs_conjugated_BIE','plus_vs_native_BIE', ...
        'minus_screen_vs_native_BIE'};
    for p=1:3
        m=metrics(fields{p}(mask),refs{p}(mask),w);
        m.case_name=c.name; m.comparison=labels{p};
        rows(3*(j-1)+p)=m;
    end
    native_m=metrics(minus(mask),native(mask),w);
    conjugated_m=metrics(conj(minus(mask)),conj(native(mask)),w);
    assert(abs(native_m.complex_l2-conjugated_m.complex_l2)<1e-14);
    % Paraxial, short-distance amplitude transport: |G_plus|-1 ~ -delta.
    % All M99 points are inside the untapered |x|<42 m surface plateau.
    x=d.x_m(:); assert(all(abs(x(mask))<42));
    eta_prime=c.amplitude_m*.1*cos(.1*x);
    eta_second=-c.amplitude_m*.1^2*sin(.1*x);
    delta=z*(eta_second+2*(amplitude_prime./amplitude).*eta_prime);
    actual=abs(plus(mask))-1;
    transport_residual=sqrt(sum(w.*(actual+delta(mask)).^2)/sum(w.*actual.^2));
    cases(j).name=c.name; cases(j).G_plus=plus; cases(j).G_minus=minus;
    cases(j).G_BIE_native=native; cases(j).mask=mask;
    cases(j).saved_plus_max_abs_difference=saved_closure;
    cases(j).transport_delta=delta;
    cases(j).transport_relative_residual=transport_residual;
    cases(j).mirror_RMS=sqrt(sum(w.*(abs(plus(mask))+abs(native(mask))-2).^2));
end
% Exact normal plane wave, physical surface z=h, time exp(-i*omega*t).
% Path change is -2h, so the reflected rough/flat ratio is exp(-2i*k*h).
h=.01; exact=exp(-2i*k*h);
plane=struct('height_m',h,'exact_ratio',exact,'plus_ratio',exp(2i*k*h), ...
    'minus_ratio',exp(-2i*k*h),'plus_absolute_error',abs(exp(2i*k*h)-exact), ...
    'minus_absolute_error',abs(exp(-2i*k*h)-exact));
assert(plane.minus_absolute_error==0 && plane.plus_absolute_error>.1);
table_rows=struct2table(rows);
audit=struct('status','HEIGHT_SIGN_MISMATCH_IN_FIXED_COMPARISON_CONFIRMED', ...
    'scope','validation-only saved-start +/- screen evaluation; not production remediation', ...
    'run_dir',out,'input_files',{{line_file,input_file}},'x_m',d.x_m, ...
    'weights',d.weights,'frequency_hz',4000,'sound_speed_mps',1500, ...
    'receiver_depth_m',z,'time_convention','exp(-i*omega*t)', ...
    'physical_surface','z=eta, z positive downward', ...
    'flat_reproduction_relative_l2',flat_closure,'plane_control',plane, ...
    'metrics',table_rows,'cases',cases,'matlab_version',version, ...
    'new_saved_start_FFT_fields',6,'new_source_propagations',0, ...
    'new_BIE_solves',0,'new_BIE_field_evaluations',0,'new_Bellhop_jobs',0);
mkdir(out);
save(fullfile(out,'height_sign_audit.mat'),'audit','-v7.3');
writetable(table_rows,fullfile(out,'metrics.csv'));
summary=rmfield(audit,{'cases','metrics','x_m','weights','plane_control'});
summary.plane_plus_absolute_error=plane.plus_absolute_error;
summary.plane_minus_absolute_error=plane.minus_absolute_error;
summary.max_saved_plus_reproduction=max([cases.saved_plus_max_abs_difference]);
summary.transport_relative_residual=[cases.transport_relative_residual];
fid=fopen(fullfile(out,'audit_summary.json'),'w','n','UTF-8');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));
disp(table_rows); disp(summary);
end

function m=metrics(a,b,w)
aa=abs(a); bb=abs(b); ac=aa-sum(w.*aa); bc=bb-sum(w.*bb);
m=struct('complex_l2',sqrt(sum(w.*abs(a-b).^2)/sum(w.*bb.^2)), ...
    'magnitude_relative_l2',sqrt(sum(w.*(aa-bb).^2)/sum(w.*bb.^2)), ...
    'phase_rms_rad',sqrt(sum(w.*angle(a.*conj(b)).^2)), ...
    'magnitude_signed_correlation',sum(w.*ac.*bc)/sqrt(sum(w.*ac.^2)*sum(w.*bc.^2)), ...
    'sample_count',numel(a));
end
