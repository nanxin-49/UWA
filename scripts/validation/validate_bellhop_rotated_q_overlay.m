function receipt = validate_bellhop_rotated_q_overlay(options)
%VALIDATE_BELLHOP_ROTATED_Q_OVERLAY Compare generated overlay with official Q.
arguments
    options.official_exe (1,:) char
    options.overlay_exe (1,:) char
    options.evidence_root (1,:) char
    options.output_dir (1,:) char
end
assert(~isfolder(options.output_dir)); mkdir(options.output_dir);
diary(fullfile(options.output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
t=tic; status='FAILED'; err=''; fields={}; processes={};
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
try
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM preflight gate.');
% New C-mode qualification: uniform analytic normalization and variable Q.
for j=1:2
    cfg=make_bellhop_rotated_munk_config('up');
    cfg.run_type='C'; cfg.ssp_spacing_m=2.5; cfg.step_m=.5; cfg.uniform_control=j==1;
    pressures=cell(1,2);
    for ib=1:2
        case_dir=fullfile(options.output_dir,sprintf('field%d_binary%d',j,ib)); mkdir(case_dir);
        case_root=fullfile(case_dir,'case');
        meta=write_bellhop_rotated_munk_q(case_root,cfg);
        if ib==2, write_mode(case_root,2); end
        exe=options.official_exe; if ib==2, exe=options.overlay_exe; end
        processes{end+1}=run_bellhop_rotated_local_case(exe,case_root,20);
        pressures{ib}=read_bellhop_shd_unfolded_vertical([case_root '.shd']);
    end
    official=pressures{1}.pressure; overlay=pressures{2}.pressure;
    add(sprintf('field%d_overlay_complex_relative',j),norm(overlay-official)/norm(official),1e-4);
    r=cfg.receiver_range_m; cs=munk_profile_sound_speed(cfg.source_z_m,cfg.munk_parameters);
    if j==1
        cs=cfg.uniform_c_mps; tau=r/cs;
        native_amp=1./r; physical_amp=native_amp;
    else
        tau=munk_vertical_travel_time_reference(cfg.source_z_m-r,cfg.source_z_m,cfg.munk_parameters);
        I=arrayfun(@(rr)integral(@(v)munk_profile_sound_speed(cfg.source_z_m-v,cfg.munk_parameters), ...
            0,rr,'AbsTol',1e-8,'RelTol',1e-12),r);
        ce=munk_profile_sound_speed(cfg.source_z_m-r,cfg.munk_parameters);
        native_amp=sqrt(ce./(r.*I));
        % Independent 3D geometrical transport in a plane-stratified medium.
        % This is a diagnostic comparison, never a fitted correction.
        physical_amp=sqrt(cs*ce)./I;
    end
    expected=-native_amp.*exp(-1i*2*pi*cfg.frequency_hz*tau);
    add(sprintf('field%d_native_transport_amplitude',j),max(abs(abs(official)./native_amp-1)),2e-3);
    add(sprintf('field%d_native_transport_phase_rad',j),max(abs(angle(official.*conj(expected)))),.05);
    fields{j}=struct('cfg',cfg,'meta',meta,'official',pressures{1},'overlay',pressures{2}, ...
        'native_transport_amplitude',native_amp,'plane_stratified_3d_transport_amplitude',physical_amp, ...
        'physical_to_native_amplitude_ratio',physical_amp./native_amp,'reference_time_s',tau, ...
        'note','Native R spreading remains diagnostic for variable Q; no correction was applied.');
end
% Replay the exact already validated raw A/R inputs with the new executable.
for id=1:2
    direction={'up','down'}; direction=direction{id};
    src=fullfile(options.evidence_root,'20261114_s1_rays_01',direction);
    folder=fullfile(options.output_dir,['ray_' direction]); mkdir(folder);
    for ext={'.env','.ssp'}, copyfile(fullfile(src,['case' ext{1}]),fullfile(folder,['case' ext{1}])); end
    case_root=fullfile(folder,'case'); write_mode(case_root,2);
    processes{end+1}=run_bellhop_rotated_local_case(options.overlay_exe,case_root,20);
    r0=ray_arrays(fullfile(src,'case.ray')); r1=ray_arrays([case_root '.ray']);
    assert(isequal(size(r0),size(r1)),'Overlay/native ray record-size difference.');
    add([direction '_overlay_ray_difference_m'],max(abs(r1(:,2:3)-r0(:,2:3)),[],'all'),1e-5);
    src=fullfile(options.evidence_root,'20261114_s1_direct_02',[direction '_refine3']);
    folder=fullfile(options.output_dir,['arrival_' direction]); mkdir(folder);
    for ext={'.env','.ssp'}, copyfile(fullfile(src,['case' ext{1}]),fullfile(folder,['case' ext{1}])); end
    case_root=fullfile(folder,'case'); write_mode(case_root,2);
    processes{end+1}=run_bellhop_rotated_local_case(options.overlay_exe,case_root,20);
    a0=read_bellhop_arrivals_ascii_vertical(fullfile(src,'case.arr'));
    a1=read_bellhop_arrivals_ascii_vertical([case_root '.arr']);
    assert(height(a0)==height(a1));
    add([direction '_overlay_delay_difference_s'],max(abs(a0.delay_s-a1.delay_s)),1e-8);
end
% Verify that this first derivative cannot silently run reflected Q physics.
folder=fullfile(options.output_dir,'reject_reflected_Q'); mkdir(folder);
cfg=make_bellhop_rotated_munk_config('up'); cfg.run_type='C';
case_root=fullfile(folder,'case'); write_bellhop_rotated_munk_q(case_root,cfg);
fid=fopen([case_root '.iwcov'],'w'); assert(fid>=0);
fprintf(fid,'1\n25\n50\n0\n3\n25 0\n25 1000\n25 2000\n'); fclose(fid);
try
    run_bellhop_rotated_local_case(options.overlay_exe,case_root,20);
catch expected_error
    fid=fopen(fullfile(folder,'expected_rejection.txt'),'w'); fprintf(fid,'%s',expected_error.message); fclose(fid);
end
prt=fileread([case_root '.prt']);
add('reflected_Q_rejected',double(~contains(prt,'POC requires uniform lossless c = 1500 m/s')),0);
status='COMPLETE';
catch exception
err=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
receipt=struct('status',status,'elapsed_s',toc(t),'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'matlab_release',version,'result_dir',options.output_dir, ...
    'official_exe_sha256',black_file_hash(options.official_exe), ...
    'overlay_exe_sha256',black_file_hash(options.overlay_exe), ...
    'scope','Direct Q overlay only. Rotated point spreading, physical reflections and total field remain unqualified.');
save(fullfile(options.output_dir,'result.mat'),'checks','receipt','fields','processes');
fid=fopen(fullfile(options.output_dir,'error.log'),'w'); fprintf(fid,'%s',err); fclose(fid);
fid=fopen(fullfile(options.output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
if ~isempty(err), fprintf('%s\n',err); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function write_mode(stem,mode)
fid=fopen([stem '.iwcov'],'w'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%d\n0\n0\n0\n0\n',mode);
end
function output=ray_arrays(file)
fid=fopen(file,'r'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fgetl(fid); fscanf(fid,'%f',1); fscanf(fid,'%d',3); beams=fscanf(fid,'%d',2);
fscanf(fid,'%f',2); fgetl(fid); fgetl(fid);
data=cell(beams(1),1);
for j=1:beams(1)
    fscanf(fid,'%f',1); info=fscanf(fid,'%d',3);
    rq=fscanf(fid,'%f',[2 info(1)]).';
    data{j}=[j*ones(size(rq,1),1),rq];
end
output=vertcat(data{:});
end