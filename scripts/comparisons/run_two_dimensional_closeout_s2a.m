function result = run_two_dimensional_closeout_s2a(run_id,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_CLOSEOUT_S2A Bounded S2a source/flat/line audit.
% No BIE solve or density evaluation, Bellhop job, PM, or x-z field run.
arguments
    run_id (1,:) char
    available_physical_bytes (1,1) double {mustBeFinite,mustBePositive}
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_-]+$','once')),'Invalid run ID.');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(root,'results','comparisons','two_dimensional_channel_closeout',run_id);
assert(~isfolder(out),'Refusing to overwrite a run.');
oldroot = fullfile(root,'results','visualization');
input_file = fullfile(oldroot,'pe_bellhop_bie_xz','20261019_115214_348', ...
    'source_data','numeric_inputs.mat');
corrected_file = fullfile(oldroot,'pe_bellhop_bie_corrected', ...
    '20261021_native_reset_v2','plot_data','comparison_data.mat');
xz_file = fullfile(oldroot,'pe_bellhop_bie_corrected', ...
    '20261021_native_reset_v2','plot_data','xz_data.mat');
fk_file = fullfile(root,'results','validation','pe_bie_full_kirchhoff', ...
    'pe_bie_full_kirchhoff_validation.mat');
s0_audit_file=fullfile(root,'results','comparisons', ...
    'two_dimensional_channel_closeout','20261020_2330_s0s1', ...
    'audit','independent_audit.json');
s0_manifest_file=fullfile(root,'results','comparisons', ...
    'two_dimensional_channel_closeout','20261020_2330_s0s1', ...
    'source_manifest.csv');
e10_file=fullfile(root,'results','validation','pe_bie_height_sign_audit', ...
    '20261021_mirror_diagnostic','height_sign_audit.mat');
e11_audit_file=fullfile(oldroot,'pe_bellhop_bie_corrected', ...
    '20261021_native_reset_v2','metrics','corrected_audit.json');
source_files = {s0_audit_file,s0_manifest_file,e10_file,e11_audit_file, ...
    input_file,corrected_file,xz_file,fk_file, ...
    fullfile(root,'scripts','validation','support','run_pe_1d_surface_reflection_validation.m'), ...
    fullfile(root,'scripts','validation','support','run_pe_physical_surface_model0_validation.m'), ...
    fullfile(root,'scripts','validation','support','fk_single_receiver_incident.m'), ...
    fullfile(root,'scripts','validation','support','solve_full_kirchhoff_2d.m'), ...
    fullfile(root,'scripts','validation','audit_pe_bie_height_phase_sign.m'), ...
    fullfile(root,'scripts','validation','prepare_pe_bellhop_bie_corrected_data.m'), ...
    fullfile(root,'scripts','validation','validate_full_kirchhoff_bie_vertical.m'), ...
    fullfile(root,'scripts','validation','support','solve_helmholtz_bie_halfplane_vertical.m'), ...
    [mfilename('fullpath') '.m'], ...
    fullfile(root,'scripts','reporting','render_two_dimensional_closeout_s2a.m'), ...
    fullfile(root,'scripts','validation','support','audit_two_dimensional_closeout_s2a.m')};
assert(all(cellfun(@isfile,source_files)),'Required input/code missing.');
s0audit=jsondecode(fileread(s0_audit_file));
e11audit=jsondecode(fileread(e11_audit_file));
assert(s0audit.pass && s0audit.figures_pixel_identical);
assert(strcmp(e11audit.status,'PASS_SAVED_DATA_AND_PIXEL_REDRAW') && ...
    e11audit.redraw_pixel_differences==0);
s0manifest=readtable(s0_manifest_file,'TextType','string');
inherited=table(s0manifest.source_id,s0manifest.relative_path, ...
    s0manifest.sha256,strings(height(s0manifest),1),false(height(s0manifest),1), ...
    'VariableNames',{'source_id','relative_path','expected_sha256', ...
    'actual_sha256','hash_matches'});
for ii=1:height(inherited)
    path=fullfile(root,char(inherited.relative_path(ii)));
    assert(isfile(path));
    inherited.actual_sha256(ii)=string(local_sha256(path));
    inherited.hash_matches(ii)=strcmpi(inherited.actual_sha256(ii),inherited.expected_sha256(ii));
end
assert(all(inherited.hash_matches),'Inherited E01-E09 evidence has drifted.');
t=load(input_file,'inputs'); in=t.inputs;
t=load(corrected_file,'data'); d=t.data;
t=load(xz_file,'data'); xz=t.data;
t=load(fk_file,'validation'); histfk=t.validation;
k=2*pi*4000/1500; N=984; assert(numel(in.xfull_m)==N);
assert(isequal(d.x_m(:),in.xfull_m(:)) && isequal(logical(d.mask(:)),logical(in.original_M99(:))));
assert(nnz(d.mask)==411 && isequal(xz.z_m(:),(.5:.125:6).'));
assert(isequal(xz.x_m(:),d.x_m(d.mask)));
assert(in.config.frequency_hz==4000 && in.config.c0_mps==1500 && ...
    in.config.z_tx_m==100 && in.config.z_rx_m==3 && ...
    in.config.sigma_src_m==.3 && in.config.nx==N && ...
    in.config.xw_m==192.1875 && in.config.pe_step_m==.05);
assert(max(abs(in.pe.source_field(:)-exp(-.5*(in.xfull_m(:)/.3).^2)))==0);
assert(max(abs(d.source_coefficients(:)-fft(in.pe.source_field(:))))<1e-12);
assert(max(abs(d.source_kx(:)-in.pe.kx_rad_per_m(:)))==0);
mask=logical(d.mask(:)); xall=d.x_m(:); x=xall(mask); z=3+zeros(size(x));
w=d.weights(mask); w=w/sum(w);
src=struct('kx',in.pe.kx_rad_per_m(:), ...
    'kz',sqrt(complex(k^2-in.pe.kx_rad_per_m(:).^2,0)), ...
    'coeff',fft(in.pe.source_field(:)),'n_original',N, ...
    'x0_m',xall(1),'z_tx_m',100);
[inc0,~,~]=fk_single_receiver_incident(src,xall,zeros(N,1));
incident_relative=norm(inc0-in.pe.surface_incident_field(:)*exp(1i*k*100))/norm(inc0);
assert(incident_relative<1e-10,'Saved PE incident differs from common spectrum.');
flat_PE=in.PE_flat_envelope(:)*exp(1i*k*103);
flat_BIE=in.BIE_flat_native(:);
flat_source_relative=norm(flat_PE-flat_BIE)/norm(flat_BIE);
assert(flat_source_relative<1e-10,'PE/BIE flat source mapping failed.');
analytic_flat=-fk_single_receiver_incident(src,xall,-3+zeros(N,1));
flat_analytic_relative=norm(flat_BIE-analytic_flat)/norm(flat_BIE);
assert(flat_analytic_relative<1e-10,'BIE flat source mapping failed.');
saved_inc_rel=zeros(2,1); saved_dn_rel=zeros(2,1); saved_normal_max=zeros(2,1);
for jj=1:2
    det=histfk.details(jj).refined;
    sx=det.surface_x_m(:); sz=det.surface_z_m(:);
    A=[.01 .20]; [~,ep]=local_surface(sx,A(jj));
    normal=[-ep./hypot(1,ep),1./hypot(1,ep)];
    [ui,ux,uz]=fk_single_receiver_incident(src,sx,sz);
    dni=normal(:,1).*ux+normal(:,2).*uz;
    saved_inc_rel(jj)=norm(ui-det.surface_incident(:))/norm(ui);
    saved_dn_rel(jj)=norm(dni-det.surface_incident_dn(:))/norm(dni);
    saved_normal_max(jj)=max(abs(normal-det.surface_normal),[],'all');
end
assert(max([saved_inc_rel;saved_dn_rel;saved_normal_max])<1e-10, ...
    'Saved FullFK input/normal does not match source contract.');

% Exercise the actual validation entry: depth adapter -> screen -> march -> ratio.
base=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',192.1875, ...
    'nx',N,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
    'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
    'reflection_model','model0_normal');
hvals=[0 .01 -.01]; plane_rows=cell(3,1); plane_results=cell(3,1);
for j=1:3
    cfg=base; cfg.source_field_override=ones(N,1); cfg.eta_z_m=hvals(j)+zeros(N,1);
    st=tic; q=run_pe_physical_surface_model0_validation(cfg); sec=toc(st);
    plane_results{j}=q;
    path_flat=100+3; path_shift=(100-hvals(j))+(3-hvals(j));
    expected=-exp(1i*k*path_shift);
    ratio=q.reflected_native./plane_results{1}.reflected_native;
    ratio_expected=exp(1i*k*(path_shift-path_flat));
    assert(max(abs(q.source_field(:)-1))==0);
    assert(max(abs(q.physical_eta_z_m(:)-hvals(j)))==0);
    assert(max(abs(q.adapter_eta_up_m(:)+hvals(j)))==0);
    assert(max(abs(q.reflected_native(:)-expected))<1e-10);
    assert(max(abs(ratio-ratio_expected))<1e-10);
    plane_rows{j}=table(hvals(j),real(ratio_expected),imag(ratio_expected), ...
        max(abs(q.reflected_native(:)-expected)),max(abs(ratio-ratio_expected)),sec, ...
        'VariableNames',{'eta_z_m','expected_G_real','expected_G_imag', ...
        'raw_field_max_abs_error','G_max_abs_error','wall_seconds'});
end
plane_control=vertcat(plane_rows{:});
cfg=base; cfg.source_field_override=in.pe.source_field(:); cfg.eta_z_m=zeros(N,1);
st=tic; peflat=run_pe_physical_surface_model0_validation(cfg); peflat_seconds=toc(st);
assert(norm(peflat.surface_incident_field(:)-in.pe.surface_incident_field(:))/ ...
    norm(in.pe.surface_incident_field(:))<1e-10);
assert(norm(peflat.reflected_field(:)-in.PE_flat_envelope(:))/norm(in.PE_flat_envelope(:))<1e-10);
PE_rough=cell(3,1); BIE_rough=cell(3,1); pe_seconds=zeros(3,1);
for j=1:3
    cfg.eta_z_m=d.cases(j).eta_m(:);
    st=tic; q=run_pe_physical_surface_model0_validation(cfg); pe_seconds(j)=toc(st);
    g=q.reflected_field(:)./peflat.reflected_field(:);
    assert(max(abs(g(mask)-d.cases(j).G_PE(mask)))<1e-10, ...
        'Corrected saved PE line differs from actual adapter entry.');
    PE_rough{j}=d.cases(j).G_PE(:).*flat_PE;
    BIE_rough{j}=d.cases(j).G_BIE_native(:).*flat_BIE;
end
assert(max(abs(BIE_rough{3}(mask)-xz.BIE_native(:,xz.z_m==3)))<1e-11);
assert(max(abs(PE_rough{3}(mask)-xz.PE_carrier_restored(:,xz.z_m==3)))<1e-10);

% S2a representative FullFK cost, with 30% physical-RAM reserve.
% The caller samples Windows \Memory\Available Bytes immediately before
% invocation. MATLAB memory() may include the page file on this host.
available=available_physical_bytes;
% Conservative live-block estimate: twelve complex receiver-by-surface
% matrices plus 20 MiB graphics/FFT overhead. Source evaluation is blocked
% and finishes before the receiver Green blocks are allocated.
estimated_peak=ceil(numel(x)*4097*16*12*1.2)+20*1024^2;
assert(estimated_peak<=.7*available,'FullFK memory preflight failed.');
fkout=cell(4,1); fk_seconds=zeros(4,1); fk_reuse=cell(4,1);
fk_saved={histfk.flat_refined(:),histfk.details(1).refined.receiver_field(:), ...
    [],histfk.details(2).refined.receiver_field(:)};
fk_saved_relative=nan(4,1); fk_term_closure=nan(4,1);
amplitudes=[0 .01 .05 .20];
for j=1:4
    [eta,~]=local_surface(xall,amplitudes(j));
    if j>1
        assert(max(abs(eta-d.cases(j-1).eta_m(:)))<1e-12);
    end
    if j==1 || j==2 || j==3
        % Flat/weak are representative equality checks; medium is missing.
        st=tic; fkout{j}=local_fk(src,amplitudes(j),x,z); fk_seconds(j)=toc(st);
        assert(fk_seconds(j)<900,'Single FullFK job exceeded 900 s.');
        if ~isempty(fk_saved{j})
            fk_saved_relative(j)=norm(fkout{j}.receiver_field-fk_saved{j}(mask))/ ...
                norm(fk_saved{j}(mask));
            assert(fk_saved_relative(j)<1e-10,'Saved FullFK mismatch.');
            fk_reuse{j}='saved E12 verified by current full integral';
        else
            fk_reuse{j}='new medium line';
        end
    else
        % Strong source, geometry, normals, and code version were checked
        % against the same weak and flat saved E12 implementation.
        q=histfk.details(2).refined;
        assert(isequal(q.receiver_x_m(:),xall) && ...
            max(abs(q.surface_z_m-local_surface(q.surface_x_m,.20)))<1e-12);
        fkout{j}=struct('receiver_field',fk_saved{j}(mask), ...
            'receiver_term_incident_dn',nan(size(x)), ...
            'receiver_term_green_dn',nan(size(x)));
        fk_reuse{j}='saved E12 strong line; current flat/weak code equality';
    end
    if j<4
        fk_term_closure(j)=norm(fkout{j}.receiver_field - ...
            (fkout{j}.receiver_term_incident_dn+fkout{j}.receiver_term_green_dn));
        assert(fk_term_closure(j)<1e-10);
    end
end
assert(sum([plane_control.wall_seconds])+peflat_seconds+sum(pe_seconds)+sum(fk_seconds)<7200);

names={'flat','weak_low_K','region_II_low_K','strong_height_low_K'};
PE={flat_PE(mask),PE_rough{1}(mask),PE_rough{2}(mask),PE_rough{3}(mask)};
BIE={flat_BIE(mask),BIE_rough{1}(mask),BIE_rough{2}(mask),BIE_rough{3}(mask)};
FK=cellfun(@(a)a.receiver_field,fkout,'UniformOutput',false);
fields=struct('x_m',x,'z_m',z,'weights',w,'case_names',{names}, ...
    'PE',{PE},'BIE',{BIE},'FullFK',{FK}, ...
    'FK_components',{fkout},'FullFK_provenance',{fk_reuse}, ...
    'source',src,'plane_control',plane_control);
rows=cell(4*3,1); n=0;
for j=1:4
    pairs={{'PE','BIE'},{'FullFK','BIE'},{'PE','FullFK'}};
    for p=1:3
        n=n+1; aa=fields.(pairs{p}{1}){j}; bb=fields.(pairs{p}{2}){j};
        rows{n}=local_metric_row(names{j},pairs{p}{1},pairs{p}{2},aa,bb,w);
    end
end
metrics=vertcat(rows{:});
line_rows=cell(4*411*3,1); n=0;
for j=1:4
    for m=1:3
        methods={'PE','BIE','FullFK'}; u=fields.(methods{m}){j};
        for ii=1:numel(x)
            n=n+1; line_rows{n}=table(string(names{j}),string(methods{m}), ...
                ii,x(ii),3,real(u(ii)),imag(u(ii)),abs(u(ii)),angle(u(ii)),w(ii), ...
                string(local_quantity(j)),"GAUSS_X_4K_UNIT_PEAK_N984", ...
                "physical_map","S2a-depth-positive-native-v1", ...
                "exp(-i omega t); outgoing exp(+ikR)", ...
                "unit-source reflected response; not Pa", ...
                "BIE_native_same_source", ...
                'VariableNames',{'case_name','method','line_index','x_m','z_m', ...
                'u_real','u_imag','u_abs','u_phase_rad','weight', ...
                'quantity_kind','source_id','source_tier','source_map_version', ...
                'phase_convention','normalization','reference_id'});
        end
    end
end
line_fields=vertcat(line_rows{:});
diagnostic_rows=cell(3*3*411,1); n=0;
for j=2:4
    for m=1:3
        methods={'PE','BIE','FullFK'};
        flat=fields.(methods{m}){1}; rough=fields.(methods{m}){j};
        g=rough./flat; increment=rough-flat;
        for ii=1:numel(x)
            n=n+1;
            diagnostic_rows{n}=table(string(names{j}),string(methods{m}),ii,x(ii), ...
                real(g(ii)),imag(g(ii)),real(increment(ii)),imag(increment(ii)), ...
                "GAUSS_X_4K_UNIT_PEAK_N984","S2a-depth-positive-native-v1", ...
                'VariableNames',{'case_name','method','line_index','x_m', ...
                'G_real','G_imag','delta_u_real','delta_u_imag', ...
                'source_id','source_map_version'});
        end
    end
end
diagnostics=vertcat(diagnostic_rows{:});
increment_rows=cell(9,1); n=0;
for j=2:4
    pairs={{'PE','BIE'},{'FullFK','BIE'},{'PE','FullFK'}};
    for p=1:3
        n=n+1; ca=pairs{p}{1}; rb=pairs{p}{2};
        eflat=fields.(ca){1}-fields.(rb){1};
        erough=fields.(ca){j}-fields.(rb){j};
        edelta=(fields.(ca){j}-fields.(ca){1})- ...
            (fields.(rb){j}-fields.(rb){1});
        den=sum(w.*abs(fields.(rb){1}).^2);
        flat_term=sum(w.*abs(eflat).^2)/den;
        delta_term=sum(w.*abs(edelta).^2)/den;
        cross_term=2*real(sum(w.*eflat.*conj(edelta)))/den;
        rough_term=sum(w.*abs(erough).^2)/den;
        closure=rough_term-flat_term-delta_term-cross_term;
        assert(abs(closure)<1e-10);
        increment_rows{n}=table(string(names{j}),string(ca),string(rb), ...
            flat_term,delta_term,cross_term,rough_term,closure, ...
            'VariableNames',{'case_name','candidate','reference', ...
            'flat_squared_term','delta_squared_term','cross_term', ...
            'rough_squared_term','closure_residual'});
    end
end
increment_decomposition=vertcat(increment_rows{:});
% S2b core 108 sample cells: only actual same-coordinate saved values count.
ids=strings(9,1); xx=zeros(9,1); zz=zeros(9,1); n=0;
for iz=1:3
    for ix=1:3
        n=n+1; ids(n)=sprintf('P%d%d',iz,ix);
        xx(n)=(-.1953125)+(ix-1)*.1953125;
        zz(n)=2.875+(iz-1)*.125;
    end
end
gap_rows=cell(108,1); n=0; methods={'PE','BIE','FullFK'};
for m=1:3
    for j=1:4
        for ii=1:9
            n=n+1; on_line=zz(ii)==3;
            strong_raw=(j==4 && m<=2);
            available_point=on_line || strong_raw;
            reason="saved z=3 line";
            if strong_raw && ~on_line, reason="saved strong x-z raw field"; end
            if ~available_point, reason="off-depth value absent; S2b"; end
            gap_rows{n}=table(string(methods{m}),string(names{j}),4000,ids(ii), ...
                xx(ii),zz(ii),"u_reflected_native",available_point,reason, ...
                "GAUSS_X_4K_UNIT_PEAK_N984","S2a-depth-positive-native-v1", ...
                'VariableNames',{'method','case_name','frequency_hz','receiver_id', ...
                'x_m','z_m','quantity_kind','available','reason', ...
                'source_id','convention_version'});
        end
    end
end
s2b_coverage=vertcat(gap_rows{:});
assert(height(s2b_coverage)==108 && nnz(s2b_coverage.available)==48);
s2b_rows=cell(12,1); n=0;
for m=1:3
    for j=1:4
        n=n+1;
        covered=s2b_coverage.method==string(methods{m}) & ...
            s2b_coverage.case_name==string(names{j});
        missing=nnz(~s2b_coverage.available(covered));
        method=string(methods{m});
        cost_basis="unmeasured; budget before job";
        if m==1, cost_basis="saved-start PE depth evaluations; use measured S2a ~0.3 s line as context, not estimate"; end
        if m==3, cost_basis="FullFK 4097-source-sample line measured 0.39-0.55 s; batch six points, preflight again"; end
        s2b_rows{n}=table(method,string(names{j}),missing,cost_basis, ...
            'VariableNames',{'method','case_name','missing_complex_values','cost_basis'});
    end
end
s2b_cost_tasks=vertcat(s2b_rows{:});
source_qualification=table(["PE-BIE";"FullFK-BIE";"PE-FullFK";"Bellhop-BIE"], ...
    ["SOURCE_MATCHED";"SOURCE_MATCHED";"SOURCE_MATCHED";"SOURCE_UNRESOLVED"], ...
    ["same unit Gaussian FFT source; incident and flat closure"; ...
    "same FFT coefficients/derivative and physical normal; flat checked"; ...
    "same saved spectrum, unit peak and native carrier"; ...
    "X directivity pattern only; absolute strength/phase map absent"], ...
    ["native carrier restoration only";"native identity"; ...
    "native carrier restoration only";"unresolved X-to-Gaussian map"], ...
    ["no fitted scalar";"no fitted scalar";"no fitted scalar"; ...
    "fixed flat fit remains diagnostic, not applied"], ...
    'VariableNames',{'method_pair','qualification','evidence', ...
    'physical_mapping','calibration_tier'});
invalidations=table(["legacy PE plus-screen";"legacy one-sided conjugated G"; ...
    "corrected v2 PE line";"E12 FullFK flat/weak/strong"; ...
    "Bellhop native SHD"], ...
    ["diagnostic only";"diagnostic only";"reused after actual-entry closure"; ...
    "reused after current flat/weak equality";"source unresolved"], ...
    'VariableNames',{'artifact','S2a_disposition'});
source_manifest=table(string(source_files(:)),strings(numel(source_files),1), ...
    'VariableNames',{'path','sha256'});
for ii=1:numel(source_files)
    source_manifest.sha256(ii)=string(local_sha256(source_files{ii}));
end
cost=struct('core_limit_seconds',7200,'single_job_limit_seconds',900, ...
    'physical_memory_available_preflight_bytes',available, ...
    'physical_memory_measurement','Windows performance counter Available Bytes', ...
    'estimated_peak_increment_bytes',estimated_peak, ...
    'reserve_fraction',.30,'plane_seconds',sum(plane_control.wall_seconds), ...
    'PE_flat_seconds',peflat_seconds,'PE_rough_seconds',pe_seconds, ...
    'FullFK_seconds',fk_seconds,'new_FullFK_lines',1, ...
    'representative_FullFK_rechecks',2,'new_BIE_solves',0, ...
    'new_BIE_density_evaluations',0,'new_Bellhop_jobs',0, ...
    'new_full_xz_fields',0);
checks=struct('incident_PE_vs_spectrum_relative',incident_relative, ...
    'PE_BIE_flat_source_relative',flat_source_relative, ...
    'BIE_flat_vs_analytic_relative',flat_analytic_relative, ...
    'plane_control_max_error',max(plane_control.G_max_abs_error), ...
    'saved_FullFK_incident_relative',saved_inc_rel, ...
    'saved_FullFK_normal_derivative_relative',saved_dn_rel, ...
    'saved_FullFK_normal_max_abs',saved_normal_max, ...
    'FullFK_saved_relative',fk_saved_relative, ...
    'FullFK_component_closure',fk_term_closure, ...
    'coverage_cells_available',nnz(s2b_coverage.available), ...
    'coverage_cells_missing',nnz(~s2b_coverage.available), ...
    'convention_version','S2a-depth-positive-native-v1');
plot_data=struct('x_m',x,'weights',w,'case_names',{names}, ...
    'PE',{PE},'BIE',{BIE},'FullFK',{FK},'plane_control',plane_control);
mkdir(out); mkdir(fullfile(out,'data'));
save(fullfile(out,'data','s2a_data.mat'),'fields','metrics','checks','cost', ...
    'source_qualification','invalidations','s2b_coverage', ...
    'diagnostics','increment_decomposition','s2b_cost_tasks','inherited','-v7.3');
save(fullfile(out,'data','plot_data.mat'),'plot_data','-v7.3');
writetable(source_manifest,fullfile(out,'data','source_manifest.csv'));
writetable(inherited,fullfile(out,'data','e01_e09_source_recheck.csv'));
writetable(source_qualification,fullfile(out,'data','source_qualification.csv'));
writetable(invalidations,fullfile(out,'data','version_invalidations.csv'));
writetable(plane_control,fullfile(out,'data','plane_control.csv'));
writetable(metrics,fullfile(out,'data','line_metrics.csv'));
writetable(line_fields,fullfile(out,'data','line_fields.csv'));
writetable(diagnostics,fullfile(out,'data','G_delta_u.csv'));
writetable(increment_decomposition,fullfile(out,'data','increment_decomposition.csv'));
writetable(s2b_coverage,fullfile(out,'data','s2b_coverage_gap.csv'));
writetable(s2b_cost_tasks,fullfile(out,'data','s2b_cost_tasks.csv'));
local_json(fullfile(out,'data','cost.json'),cost);
local_json(fullfile(out,'data','checks.json'),checks);
render_two_dimensional_closeout_s2a(fullfile(out,'data','plot_data.mat'), ...
    fullfile(out,'figures'));
result=struct('run_dir',out,'checks',checks,'cost',cost, ...
    'source_qualification',source_qualification);
end

function q=local_quantity(j)
if j==1, q='u_flat_reflected'; else, q='u_rough_reflected'; end
end

function out=local_fk(src,A,xr,zr)
sx=linspace(-60,60,4097).';
[eta,ep]=local_surface(sx,A);
[ui,ux,uz]=fk_single_receiver_incident(src,sx,eta);
jac=hypot(1,ep); nv=[-ep./jac,1./jac];
dni=nv(:,1).*ux+nv(:,2).*uz;
cfg=struct('frequency_hz',4000,'c0_mps',1500, ...
    'surface_x_m',sx,'surface_z_m',eta,'surface_normal',nv, ...
    'surface_incident',ui,'surface_incident_dn',dni, ...
    'receiver_x_m',xr,'receiver_z_m',zr, ...
    'quadrature_measure','arc_length','local_radius_m',Inf);
out=solve_full_kirchhoff_2d(cfg);
assert(all(isfinite(out.receiver_field)) && ...
    strcmp(out.quadrature_measure,'arc_length') && isinf(out.local_radius_m));
end

function [eta,ep]=local_surface(x,A)
r=abs(x); chi=zeros(size(x)); dc=zeros(size(x));
chi(r<=42)=1; mid=r>42&r<50;
t=(r(mid)-42)/8;
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
dc(mid)=(-30*t.^2+60*t.^3-30*t.^4)/8;
eta=A*sin(.1*x).*chi;
ep=A*(.1*cos(.1*x).*chi+sin(.1*x).*dc.*sign(x));
end

function row=local_metric_row(case_name,candidate,reference,a,b,w)
delta=a-b; aa=abs(a); bb=abs(b); phase=angle(a.*conj(b));
den=sum(w.*bb.^2); e=sqrt(sum(w.*abs(delta).^2)/den);
mag=sum(w.*(aa-bb).^2)/den;
ph=sum(w.*4.*aa.*bb.*sin(phase/2).^2)/den;
assert(abs(e^2-mag-ph)<1e-10);
row=table(string(case_name),string(candidate),string(reference),numel(a), ...
    e,sqrt(mag),sqrt(sum(w.*phase.^2)), ...
    sqrt(sum(w.*(20*log10(aa./bb)).^2)),mag,ph,e^2-mag-ph, ...
    'VariableNames',{'case_name','candidate','reference','n_valid', ...
    'complex_l2','magnitude_l2','phase_rms_rad','db_rms', ...
    'magnitude_squared_term','phase_squared_term','decomposition_residual'});
end

function h=local_sha256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
while ~feof(fid), md.update(fread(fid,1024*1024,'*uint8')); end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function local_json(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
