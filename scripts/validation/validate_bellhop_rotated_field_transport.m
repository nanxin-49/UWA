function receipt = validate_bellhop_rotated_field_transport(options)
%VALIDATE_BELLHOP_ROTATED_FIELD_TRANSPORT Fixed small 3D/precision/search pilot.
arguments
    options.overlay_exe (1,:) char
    options.output_dir (1,:) char
    options.reuse_dir (1,:) char = ''
end
assert(~isfolder(options.output_dir)); mkdir(options.output_dir); t=tic;
diary(fullfile(options.output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
tags={'uniform_native','uniform_point','munk_linear_native','munk_fast_native','munk_point'};
point=[false true false false true]; fast=[true true false true true]; runs=cell(1,5);
status='FAILED'; err='';
try
for j=1:5
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM gate.');
    assert(toc(t)<60,'Transport pilot cumulative time gate.');
    cfg=make_bellhop_rotated_munk_config('up'); cfg.run_type='C';
    cfg.uniform_control=j<=2; cfg.ssp_spacing_m=2.5; cfg.step_m=.5; cfg.beam_count=1001;
    cfg.receiver_range_m=[50.003 1950.007 3950.009 4940.013];
    folder=fullfile(options.output_dir,tags{j}); mkdir(folder); stem=fullfile(folder,'case');
    meta=write_bellhop_rotated_munk_q(stem,cfg);
    fid=fopen([stem '.iwcov'],'w'); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
    write_bellhop_rotated_field_options(stem,cfg.receiver_range_m,meta.q_receiver_m,point3d=point(j),fast_q=fast(j));
    reuse=false;
    if ~isempty(options.reuse_dir)
        prior_cp=jsondecode(fileread(fullfile(options.reuse_dir,'checkpoint.json')));
        assert(strcmp(prior_cp.overlay_exe_sha256,black_file_hash(options.overlay_exe)));
        prior=fullfile(options.reuse_dir,tags{j});
        assert(strcmp(strtrim(fileread(fullfile(prior,'state.txt'))),'TERMINAL'));
        assert(str2double(fileread(fullfile(prior,'exit_code.txt')))==0);
        for ext={'.env','.ssp','.iwcov','.r25d'}
            assert(strcmp(black_file_hash([stem ext{1}]),black_file_hash(fullfile(prior,['case' ext{1}]))));
        end
        for name={'case.shd','case.prt','case.r25diag','case.r25sample','stdout.log','stderr.log','exit_code.txt','state.txt'}
            copyfile(fullfile(prior,name{1}),fullfile(folder,name{1}));
            assert(strcmp(black_file_hash(fullfile(prior,name{1})),black_file_hash(fullfile(folder,name{1}))));
        end
        process=struct('reused_from',prior,'status','TERMINAL_WITH_MATCHING_INPUT_AND_BINARY_SHA');
        reuse=true;
    end
    if ~reuse
        process=run_bellhop_rotated_local_case(options.overlay_exe,stem,min(20,60-toc(t)));
    end
    field=read_bellhop_rotated_field(stem);
    d=readmatrix([stem '.r25diag'],'FileType','text','CommentStyle','#');
    d=d(~all(isnan(d),2),:); assert(size(d,2)==10);
    pq=sind(d(:,1))./d(:,7);
    error=(d(:,3)-cfg.transverse_origin_m)-pq.*d(:,6);
    add([tags{j} '_outplane_invariant'],max(abs(error))/max(1,max(abs(d(:,3)-cfg.transverse_origin_m))),1e-9);
    [~,ia]=min(abs(d(:,1)));
    R=d(ia,2); rr=meta.range_m; cc=meta.c_mps; nodes=[rr(rr<R);R];
    expected_I=trapz(nodes,interp1(rr,cc,nodes));
    add([tags{j} '_axis_integral_relative'],abs(d(ia,6)-expected_I)/expected_I,1e-10);
    add([tags{j} '_ray_unit_error'],max(d(:,8)),1e-7);
    runs{j}=struct('cfg',cfg,'meta',meta,'field',field,'process',process,'ray_diagnostics',d);
    save(fullfile(folder,'result.mat'),'cfg','meta','field','process','d');
end
% Fixed normalization/phasor bridge; not derived by fitting these data.
r=runs{1}.cfg.receiver_range_m; frequency=runs{1}.cfg.frequency_hz;
uniform=exp(1i*2*pi*frequency*r/runs{1}.cfg.uniform_c_mps)./r;
for j=1:2
    pressure=-conj(runs{j}.field.pressure);
    add([tags{j} '_uniform_Green_complex'],norm(pressure-uniform)/norm(uniform),1e-6);
    add([tags{j} '_uniform_phase_rad'],max(abs(angle(pressure.*conj(uniform)))),1e-6);
end
add('uniform_point_correction_identity',norm(runs{1}.field.pressure-runs{2}.field.pressure)/norm(runs{1}.field.pressure),1e-6);
add('fast_Q_complex_equivalence',norm(runs{3}.field.pressure-runs{4}.field.pressure)/norm(runs{3}.field.pressure),1e-6);
cfg=runs{5}.cfg; pp=cfg.munk_parameters; cs=munk_profile_sound_speed(cfg.source_z_m,pp);
I=arrayfun(@(R)integral(@(v)munk_profile_sound_speed(cfg.source_z_m-v,pp),0,R, ...
    'RelTol',1e-12,'AbsTol',1e-8),r);
ce=munk_profile_sound_speed(cfg.source_z_m-r,pp);
amplitude=sqrt(cs*ce)./I; corrected=-conj(runs{5}.field.pressure);
add('Munk_3D_transport_amplitude',max(abs(abs(corrected)./amplitude-1)),2e-6);
tau=munk_vertical_travel_time_reference(cfg.source_z_m-r,cfg.source_z_m,pp);
add('Munk_3D_transport_phase_rad',max(abs(angle(corrected.*exp(-1i*2*pi*frequency*tau)))),.05);
ratio=abs(runs{5}.field.pressure./runs{4}.field.pressure);
samples=readmatrix(fullfile(options.output_dir,'munk_point','case.r25sample'), ...
    'FileType','text','CommentStyle','#');
samples=samples(~all(isnan(samples),2),:);
assert(size(samples,1)==numel(r) && size(samples,2)==11);
[~,order]=sort(samples(:,2)); samples=samples(order,:);
derived_factor=sqrt(samples(:,9).*samples(:,2)./(abs(cosd(samples(:,1))).*samples(:,8)));
expected_ratio=(derived_factor.*sqrt(samples(:,11)./samples(:,6))).';
add('logged_point_factor_replay',max(abs(derived_factor-samples(:,10))),1e-12);
add('logged_ray_sound_speed_replay_mps', ...
    max(abs(samples(:,5)+samples(:,4).*(samples(:,6)-samples(:,5))-samples(:,11))),1e-10);
add('Munk_point_ratio_and_receiver_c_derived',max(abs(ratio-expected_ratio)),2e-6);
metadata_error=max(abs(runs{1}.field.receiver_range_m-runs{1}.field.receiver_range_metadata_m));
add('nonzero_float32_metadata_quantization',double(metadata_error==0),0);
status='COMPLETE';
catch exception
err=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'elapsed_s',toc(t),'matlab_release',version,'result_dir',options.output_dir, ...
    'scope','Fixed direct-axis 3D transport, exact-coordinate and equivalent-search pilot; no complete field', ...
    'overlay_exe_sha256',black_file_hash(options.overlay_exe));
save(fullfile(options.output_dir,'result.mat'),'receipt','checks','runs');
fid=fopen(fullfile(options.output_dir,'error.log'),'w'); fprintf(fid,'%s',err); fclose(fid);
fid=fopen(fullfile(options.output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
if ~isempty(err), fprintf('%s\n',err); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end