function execute_p1_ds4062(project_root, run_dir)
% One formal C job per path, with a fail-closed gate before the next job.
addpath(fullfile(project_root,'scripts','validation','support'));
addpath(fullfile(run_dir,'code'),'-begin');
assert(license('test','MATLAB')==1,'MATLAB license gate');
fprintf('MATLAB_LICENSE_PASS %s\n',version);
f=4062.5;
paths={'D','S'};
original=load(fullfile(run_dir,'original_plan.mat'),'slots','x','v');
assert(isequal(original.x,[0 10 50 100 200]));
x=original.x(:);
lengths=[1484 1496];
for j=1:2
    label=paths{j};
    slot=original.slots(2*j);
    cfg=slot.cfg;
    assert(strcmp(slot.path,label) && slot.frequency_hz==f && ...
        strcmp(slot.mode,'C') && slot.beams==201 && cfg.frequency_hz==f);
    assert(cfg.source_z_m==1490 && cfg.water_depth_m==1500 && ...
        cfg.receiver_range_m==lengths(j) && cfg.step_m==.5 && cfg.ssp_spacing_m==.625);
    assert(isequal(cfg.angle_limits_deg,[-1 1]*1.1*atan2d(200,lengths(j))));
    assert(strcmpi(slot.binary_sha256,'5d87293a27e7fafc3aedefb569f3f35065294f42426503c4881d6eb55e6769e9'));
    slots=slot;
    slots.cfg.beam_count=801;
    slots.beams=801;
    slotDirs{j}=fullfile(run_dir,label); %#ok<AGROW>
    slots.output_dir=fullfile(slotDirs{j},'case_data');
    slots.local_dir=slots.output_dir;
    slots.binary=fullfile(project_root,'results','validation','bellhop_rotated_fig12b', ...
        '20261116_linux_C_01','bin','bellhop_rotated_receiver64_2020');
    save(fullfile(slotDirs{j},'plan.mat'),'slots');
    roots{j}=fullfile(slots.output_dir,'case'); %#ok<AGROW>
    order{j}=cfg.receiver_x_order; %#ok<AGROW>
    assert(isequal(sort(cfg.receiver_x_order(:)).',(1:5)));
end
baseline=load(fullfile(run_dir,'baseline4000_comparison.mat'), ...
    'x','pD','pS','pDS','refDS','EDS');
assert(isequal(baseline.x(:),x) && isequal(baseline.pDS,baseline.pD+baseline.pS));
runs=cell(1,2);
p=complex(zeros(5,2));
for j=1:2
    if j==1
        runs{j}=run_bellhop_1500m_d801_case(project_root,slotDirs{j},160);
    else
        runs{j}=run_bellhop_1500m_s801_case(project_root,slotDirs{j},160);
    end
    r=runs{j};
    assert(r.started_slots==1 && r.completed_slots==1 && r.failed_checks==0, ...
        '%s slot failed. See checkpoint_C.json.',paths{j});
    component=load(fullfile(fileparts(roots{j}),'component.mat'),'p','proc');
    raw=read_bellhop_rotated_field(roots{j});
    native=complex(zeros(5,1));
    fid=fopen([roots{j} '.shd'],'rb');
    assert(fid>=0);
    closeSHD=onCleanup(@() fclose(fid));
    for k=1:5
        assert(fseek(fid,(raw.range_record+k)*raw.record_bytes,'bof')==0);
        z=fread(fid,2,'float32');
        native(k)=complex(z(1),z(2));
    end
    clear closeSHD
    p(order{j},j)=-conj(native);
    assert(isequal(native,raw.pressure(:)) && isequal(p(:,j),component.p(:)), ...
        'Independent SHD and receiver permutation audit');
end
geom=audit_surface_event(roots{2},run_dir);
assert(all(geom.pass),'Frozen independent surface geometry gate failed');
LD=hypot(x,1484);
LS=hypot(x,1496);
refD=exp(1i*2*pi*f/1500*LD)./LD;
refS=-exp(1i*2*pi*f/1500*LS)./LS;
pD=p(:,1);
pS=p(:,2);
pDS=pD+pS;
refDS=refD+refS;
ES=abs(refS);
EDS=abs(refD)+ES;
P=[pD pS pDS];
R=[refD refS refDS];
E=[abs(refD) ES EDS];
labels=["D";"S";"D_plus_S"];
support=abs(P)>=.01*E & abs(R)>=.01*E;
residual=abs(P-R)./E;
phase=angle(P.*conj(R));
maskedPhase=abs(phase);
maskedPhase(~support)=NaN;
l2=(vecnorm(P-R)./vecnorm(R)).';
maxDE=max(residual,[],1).';
maxPhase=max(maskedPhase,[],1,'omitnan').';
metric=table(labels,l2,maxDE,maxPhase,sum(support,1).', ...
    l2<=.001,maxDE<=.002,maxPhase<=.01, ...
    'VariableNames',{'field','complex_L2','max_complex_over_E','max_phase_rad', ...
    'phase_support','L2_pass','max_complex_pass','phase_pass'});
writetable(metric,fullfile(run_dir,'DS4062_metrics.csv'));
fields=table(repmat(x,3,1),repelem(labels,5),real(P(:)),imag(P(:)), ...
    real(R(:)),imag(R(:)),E(:),residual(:),phase(:),support(:), ...
    'VariableNames',{'x_m','field','p_real','p_imag','reference_real','reference_imag', ...
    'independent_envelope','complex_error_over_E','raw_phase_error_rad','phase_support'});
csvPath=fullfile(run_dir,'DS4062_fields.csv');
writetable(fields,csvPath);
reload=readtable(csvPath,'TextType','string');
assert(isequal(reload.x_m,fields.x_m) && isequal(reload.field,fields.field) && ...
    max(abs(complex(reload.p_real,reload.p_imag)-P(:)))<1e-17,'Saved field CSV audit');
axisCheck=table([4000;f],[abs(baseline.pDS(1))/baseline.EDS(1);abs(pDS(1))/EDS(1)], ...
    [abs(baseline.refDS(1))/baseline.EDS(1);abs(refDS(1))/EDS(1)], ...
    [angle(baseline.pS(1)*conj(baseline.pD(1)));angle(pS(1)*conj(pD(1)))], ...
    'VariableNames',{'frequency_hz','axis_numeric_magnitude_over_E', ...
    'axis_analytic_magnitude_over_E','axis_S_minus_D_phase_rad'});
writetable(axisCheck,fullfile(run_dir,'axis_interference.csv'));
writetable(table(x,1490*x/1496,LD/1500,LS/1500, ...
    'VariableNames',{'x_m','analytic_reflection_x_m','D_delay_s','S_delay_s'}), ...
    fullfile(run_dir,'five_link_analytic_geometry.csv'));
save(fullfile(run_dir,'comparison.mat'),'x','LD','LS','pD','pS','pDS','refD','refS', ...
    'refDS','ES','EDS','metric','fields','geom','axisCheck','runs','baseline');
assert(all(metric.L2_pass & metric.max_complex_pass & metric.phase_pass), ...
    'Frozen field gate failed');
data=struct('x',x,'pS',pS,'refS',refS,'ES',ES,'pDS',pDS,'refDS',refDS,'EDS',EDS);
audit=plot_p1_ds_control(data,fullfile(run_dir,'S4062_DS_control'), ...
    'U0 at 4062.5 Hz | D801 + S801 | pressure-release control');
receipt=struct('status','D_S_DS_NUMERICAL_PASS_FIG_PASS', ...
    'native_jobs',2,'native_exit_codes',[runs{1}.results.process.exit_code runs{2}.results.process.exit_code], ...
    'native_elapsed_s',[runs{1}.results.process.elapsed_s runs{2}.results.process.elapsed_s], ...
    'figure',audit,'failed_gates',0,'matlab_release',version,'axis',table2struct(axisCheck), ...
    'surface_checks',height(geom),'max_R_plus1',geom.value(geom.check=="pressure_release_R_minus1"), ...
    'max_reflection_x_error_m',geom.value(geom.check=="reflection_point_x_m"));
fid=fopen(fullfile(run_dir,'comparison_receipt.json'),'w','n','UTF-8');
assert(fid>=0);
closeReceipt=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(metric);
disp(axisCheck);
end
