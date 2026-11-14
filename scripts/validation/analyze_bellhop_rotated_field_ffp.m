function receipt = analyze_bellhop_rotated_field_ffp(field_file,ffp_file,output_dir,options)
%ANALYZE_BELLHOP_ROTATED_FIELD_FFP Saved exact-grid complex comparison only.
% No solver, interpolation, fitted source phase/gain, or TL summation.
arguments
    field_file (1,:) char
    ffp_file (1,:) char
    output_dir (1,:) char
    options.scope (1,:) char {mustBeMember(options.scope,{'pilot','full'})}='pilot'
    options.suite_checkpoint (1,:) char=''
    options.ffp_sha256 (1,:) char='8e123af0beff9904289a5419fdb82f74473beda075792b7b2f9819bd133ce3be'
end
assert(~isfolder(output_dir));
[~,mem]=memory; needed=2^30; if strcmp(options.scope,'full'), needed=1.5*2^30; end
assert(mem.PhysicalMemory.Available>=needed,'Saved-field analysis RAM gate.');
assert(strcmpi(black_file_hash(ffp_file),options.ffp_sha256),'Independent FFP SHA mismatch.');
a=load(field_file,'total','envelope','plan');
assert(isfield(a,'envelope'),'Use an independently reconstructed field with component envelope.');
plan=a.plan; z=plan.depths_m(:); x=plan.x_m(:).';
if strcmp(options.scope,'pilot')
    assert(numel(z)*numel(x)<=10000,'Pilot cannot certify a complete grid.');
else
    assert(isfile(options.suite_checkpoint),'Full comparison requires completed numerical suite evidence.');
    suite=jsondecode(fileread(options.suite_checkpoint));
    assert(strcmp(suite.status,'NUMERICAL_SUITE_COMPLETE_FFP_AND_FIGURES_PENDING') && suite.failed_checks==0);
    last=suite.results; if iscell(last), last=last{end}; else, last=last(end); end
    assert(strcmpi(black_file_hash(field_file),black_file_hash(fullfile(last.result_dir,'result.mat'))), ...
        'Full comparison must use the qualified final control.');
end
mkdir(output_dir); t=tic; diary(fullfile(output_dir,'diary.log')); cleaner=onCleanup(@()diary('off'));
b=load(ffp_file,'output'); r=b.output;
assert(r.cfg.frequency_hz==plan.physical_contract.frequency_hz && r.cfg.source_depth_m==plan.physical_contract.source_z_m);
assert(r.cfg.water_depth_m==plan.physical_contract.water_depth_m && isequal(r.cfg.munk_parameters,plan.physical_contract.munk_parameters));
for name={'bottom_cp_mps','bottom_rho_gcm3','bottom_atten_db_wavelength'}, assert(isequal(r.cfg.(name{1}),plan.physical_contract.(name{1}))); end
assert(r.cfg.water_rho_gcm3==1,'Rotated transport assumes constant water density1.');
if strcmp(options.scope,'full')
    assert(isequal(z,r.depth_m(:)) && isequal(x,r.range_m(:).'),'Full comparison coordinates must be exactly identical.');
    iz=1:numel(z); ix=1:numel(x);
else
    [yes,iz]=ismembertol(z,r.depth_m(:),1e-9,'DataScale',1); assert(all(yes));
    [yes,ix]=ismembertol(x,r.range_m(:).',1e-9,'DataScale',1); assert(all(yes));
end
coordinate_error=max([abs(z-r.depth_m(iz));abs(x(:)-reshape(r.range_m(ix),[],1))]);
p_bellhop=a.total; envelope=a.envelope; p_ffp=r.pressure_1m(iz,ix);
ffp_green=r.green(iz,ix); reference_scale_error=max(abs(p_ffp-4*pi*ffp_green),[],'all','omitnan');
source_mask=plan.source_singularity_mask;
assert(isequal(size(p_bellhop),[numel(z),numel(x)]) && isequal(isnan(p_bellhop),source_mask));
assert(isequal(isnan(p_ffp),source_mask) && all(envelope(~source_mask)>0));
valid=~source_mask; phase_valid=valid & abs(p_bellhop)>=.05*envelope & abs(p_ffp)>=.05*envelope;
difference=p_bellhop-p_ffp;
envelope_error=abs(difference)./envelope;
phase_difference=angle(p_bellhop.*conj(p_ffp)); phase_difference(~phase_valid)=NaN;
wz=weights(z); wx=weights(x(:)); cell_weights=wz*wx.';
Z=z+zeros(1,numel(x)); X=zeros(numel(z),1)+x;
distance=hypot(X,Z-plan.physical_contract.source_z_m);
regions={valid,valid & Z<=4900,valid & Z>=4900 & distance>=1,valid & X==0 & distance>=1,valid & distance<1};
names={'all','bulk_z0_4900','near_bottom_r_ge1','axis_r_ge1','source_vicinity_r_lt1'};
metrics=table(strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),false(0,1), ...
    'VariableNames',{'region','points','complex_relative_L2','cell_weighted_complex_relative', ...
    'magnitude_relative_L2','max_envelope_error','phase_points','max_phase_rad','sampled_criterion_met'});
for j=1:numel(regions)
    use=regions{j}; if ~any(use,'all'), continue; end
    d=difference(use); ref=p_ffp(use); w=cell_weights(use);
    relative=norm(d)/norm(ref); weighted=sqrt(sum(w.*abs(d).^2)/sum(w.*abs(ref).^2));
    magnitude=norm(abs(p_bellhop(use))-abs(ref))/norm(ref);
    maximum=max(envelope_error(use)); ph=phase_difference(use & phase_valid);
    maxphase=max(abs(ph)); if isempty(ph), maxphase=NaN; end
    met=relative<=.02 && weighted<=.02 && maximum<=.02 && (isempty(ph)||maxphase<=.05);
    metrics(end+1,:)={names{j},nnz(use),relative,weighted,magnitude,maximum,numel(ph),maxphase,met}; %#ok<AGROW>
end
boundary=struct('top_status','NO_SAVED_TOP_POINTS','bottom_status','NOT_EVALUABLE_SPARSE_POINTS', ...
    'bottom_operator_scope','Same saved-grid one-sided derivative comparison; not full fluid DtN closure on finite x aperture');
it=find(z==0,1); ib=find(z==r.cfg.water_depth_m,1);
if ~isempty(it)
    boundary.top_status='SAMPLED_PRESSURE_RELEASE_RESIDUAL';
    boundary.top_bellhop_max_envelope=max(abs(p_bellhop(it,:))./envelope(it,:));
    boundary.top_ffp_max_pressure=max(abs(p_ffp(it,:)));
    boundary.top_tolerance=.001;
    boundary.top_sampled_criterion_met=boundary.top_bellhop_max_envelope<=.001;
end
bottom=table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'h_m','derivative_relative_L2','max_derivative_error_over_kwE', ...
    'normal_admittance_residual_bellhop','normal_admittance_residual_ffp','max_coordinate_difference_m'});
for h=[.01 .02]
    requested=r.cfg.water_depth_m-(0:4)*h;
    [present,rows]=ismembertol(requested,z,1e-9,'DataScale',1);
    if ~all(present) || isempty(ib), continue; end
    coeff=[25 -48 36 -16 3]/(12*h);
    db=coeff*p_bellhop(rows,:); df=coeff*p_ffp(rows,:);
    cw=munk_profile_sound_speed(r.cfg.water_depth_m,r.cfg.munk_parameters); omega=2*pi*r.cfg.frequency_hz;
    alpha=r.cfg.bottom_atten_db_wavelength*log(10)/20*r.cfg.frequency_hz/r.cfg.bottom_cp_mps;
    kb=omega/r.cfg.bottom_cp_mps+1i*alpha; Y=1i*kb/r.cfg.bottom_rho_gcm3;
    scale=(omega/cw)*envelope(ib,:);
    bottom(end+1,:)={h,norm(db-df)/norm(df),max(abs(db-df)./scale), ...
        norm(db-Y*p_bellhop(ib,:))/norm(scale),norm(df-Y*p_ffp(ib,:))/norm(scale), ...
        max(abs(requested-z(rows).'))}; %#ok<AGROW>
end
if height(bottom)>0, boundary.bottom_status='SAVED_GRID_DERIVATIVE_DIAGNOSTIC'; end
near=z(z>=4900 & z<=5000); near_step=inf; if numel(near)>1, near_step=max(diff(near)); end
resolution=struct('near_bottom_depth_count',numel(near),'maximum_spacing_m',near_step, ...
    'point01m_requirement_met',numel(near)>=10001 && near_step<=.010000001, ...
    'minimum_water_wavelength_m',min(r.cfg.layer_c_mps)/r.cfg.frequency_hz, ...
    'minimum_samples_per_two_way_period',min(r.cfg.layer_c_mps)/(2*r.cfg.frequency_hz*near_step), ...
    'full_depth_note','Global 5m samples are overview only and do not resolve all wavelength-scale fringes');
axis_col=find(x==0,1);
if ~isempty(axis_col)
    axis_data=table(z,real(p_bellhop(:,axis_col)),imag(p_bellhop(:,axis_col)),real(p_ffp(:,axis_col)), ...
        imag(p_ffp(:,axis_col)),envelope(:,axis_col),phase_difference(:,axis_col),envelope_error(:,axis_col), ...
        'VariableNames',{'z_m','bellhop_real','bellhop_imag','ffp_real','ffp_imag','component_envelope','phase_difference_rad','envelope_error'});
    writetable(axis_data,fullfile(output_dir,'axis_complex.csv'));
end
scope='SPARSE_DIAGNOSTIC_NOT_COMPLETE_FIELD'; if strcmp(options.scope,'full'), scope='FULL_GRID_REFERENCE_COMPARISON_NOT_AUTOMATIC_GO'; end
reference_cfg=r.cfg; clear b r a ffp_green Z X difference cell_weights
receipt=struct('status',scope,'scope',options.scope,'elapsed_s',toc(t),'result_dir',output_dir, ...
    'numerical_suite_checkpoint',options.suite_checkpoint, ...
    'analysis_code_sha256',black_file_hash([mfilename('fullpath') '.m']), ...
    'field_file',field_file,'field_sha256',black_file_hash(field_file),'ffp_file',ffp_file, ...
    'ffp_sha256',options.ffp_sha256,'coordinate_max_difference_m',coordinate_error, ...
    'FFP_4pi_Green_max_error',reference_scale_error,'reference_failed_regions',nnz(~metrics.sampled_criterion_met), ...
    'region_count',height(metrics),'reflection_order',plan.settings.max_order,'boundary',boundary,'resolution',resolution, ...
    'criteria_note','Engineering comparison at saved samples: complex L2 and max delta/E .02; phase .05rad where both |p|>=.05E. Not paper-derived tolerances or physical model certification', ...
    'weights_note','Unweighted receiving-point norm and trapezoidal cell-weighted sampled norm both saved; neither resolves uncomputed inter-sample oscillations');
writetable(metrics,fullfile(output_dir,'reference_metrics.csv')); writetable(bottom,fullfile(output_dir,'bottom_derivative_diagnostics.csv'));
save(fullfile(output_dir,'comparison.mat'),'receipt','metrics','bottom','p_bellhop','p_ffp','envelope', ...
    'envelope_error','phase_difference','phase_valid','source_mask','z','x','plan','reference_cfg','-v7.3');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); disp(jsonencode(receipt));
end
function w=weights(v)
assert(numel(v)>1 && all(diff(v)>0)); d=diff(v); w=[d(1);d(1:end-1)+d(2:end);d(end)]/2;
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end