function receipt = validate_bellhop_rotated_munk_wall(options)
%VALIDATE_BELLHOP_ROTATED_MUNK_WALL Tiny single-wall Q geometry/state gate.
% Reflected-only field is compared with unfolded Q propagation at exact
% receiver coordinates, using an independent angle-dependent fluid coefficient.
arguments
    options.overlay_exe (1,:) char
    options.output_dir (1,:) char
    options.reuse_dir (1,:) char = ''
    options.max_elapsed_s (1,1) double {mustBePositive} = 40
end
assert(~isfolder(options.output_dir)); mkdir(options.output_dir);
diary(fullfile(options.output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
t=tic; status='FAILED'; err=''; results=cell(1,2);
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
try
for id=1:2
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM gate.');
    directions={'up','down'}; cfg=make_bellhop_rotated_munk_config(directions{id});
    cfg.run_type='C'; cfg.ssp_spacing_m=5; cfg.step_m=2; cfg.beam_count=301;
    if id==1
        cfg.reflection_wall='surface'; cfg.wall_range_m=4950; cfg.range_limit_m=9950;
        cfg.receiver_range_m=[4950 4960 5950 7950 9850]; cfg.receiver_x_m=[1 0 -1];
        cfg.angle_limits_deg=[-.2 .2];
    else
        cfg.reflection_wall='bottom'; cfg.wall_range_m=50; cfg.range_limit_m=5050;
        cfg.receiver_range_m=[50 60 150 2050 4050]; cfg.receiver_x_m=[-1 0 1];
        cfg.angle_limits_deg=[-2 2];
    end
    tag=cfg.reflection_wall;
    fields=cell(1,2); metas=cell(1,2); processes=cell(1,2);
    for ib=1:2
        assert(toc(t)<options.max_elapsed_s,'Cumulative wall pilot time gate.');
        folder=fullfile(options.output_dir,sprintf('%s_mode%d',tag,ib)); mkdir(folder);
        stem=fullfile(folder,'case'); metas{ib}=write_bellhop_rotated_munk_wall_q(stem,cfg);
        if ib==2
            % Same folded medium, no reflector, no branch isolation: independent
            % geometric continuation reference, not a second physical solver.
            fid=fopen([stem '.iwcov'],'w'); assert(fid>=0); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
        end
        reuse=false;
        if ~isempty(options.reuse_dir)
            prior_cp=jsondecode(fileread(fullfile(options.reuse_dir,'checkpoint.json')));
            assert(strcmp(prior_cp.overlay_exe_sha256,black_file_hash(options.overlay_exe)), ...
                'Cannot reuse outputs from a different binary.');
            prior=fullfile(options.reuse_dir,sprintf('%s_mode%d',tag,ib));
            if isfile(fullfile(prior,'case.shd')) && isfile(fullfile(prior,'state.txt'))
                assert(strcmp(strtrim(fileread(fullfile(prior,'state.txt'))),'TERMINAL'));
                assert(str2double(fileread(fullfile(prior,'exit_code.txt')))==0);
                for ext={'.env','.ssp','.iwcov','.rotq'}
                    assert(strcmp(black_file_hash([stem ext{1}]),black_file_hash(fullfile(prior,['case' ext{1}]))), ...
                        'Refusing to reuse changed inputs.');
                end
                for retained={'case.shd','case.prt','case.iwdiag','stdout.log','stderr.log','exit_code.txt','state.txt'}
                    copyfile(fullfile(prior,retained{1}),fullfile(folder,retained{1}));
                    assert(strcmp(black_file_hash(fullfile(prior,retained{1})),black_file_hash(fullfile(folder,retained{1}))));
                end
                processes{ib}=struct('reused_from',prior,'status','TERMINAL_WITH_MATCHING_INPUT_AND_BINARY_SHA');
                reuse=true;
            end
        end
        if ~reuse
            processes{ib}=run_bellhop_rotated_local_case(options.overlay_exe,stem,min(20,options.max_elapsed_s-toc(t)));
        end
        fields{ib}=read_bellhop_shd_unfolded_vertical([stem '.shd']);
    end
    stem=fullfile(options.output_dir,sprintf('%s_mode1',tag),'case');
    fid=fopen([stem '.iwdiag'],'r'); header=strtrim(erase(fgetl(fid),'#')); fclose(fid);
    names=strsplit(header); matrix=readmatrix([stem '.iwdiag'],'FileType','text','CommentStyle','#');
    matrix=matrix(~all(isnan(matrix),2),:);
    diag=array2table(matrix,'VariableNames',names);
    assert(height(diag)==cfg.beam_count,'Every pilot beam must cross the single physical wall.');
    add([tag '_wall_residual_m'],max(abs(diag.wall_residual)),1e-7);
    add([tag '_specular_error'],max(diag.specular_error),1e-12);
    add([tag '_proper_rotation_error'],max(diag.rotation_error),1e-12);
    add([tag '_p_rotation_error'],max(diag.p_rot_error),1e-12);
    add([tag '_q_rotation_error'],max(diag.q_rot_error),1e-12);
    % Q's exact incident-side segment derivative; no right-side substitution.
    meta=metas{1}; rr=meta.range_m; cc=meta.c_mps;
    iw=find(rr==cfg.wall_range_m,1); assert(~isempty(iw) && iw>1);
    grad=(cc(iw)-cc(iw-1))/(rr(iw)-rr(iw-1)); cw=cc(iw);
    expected_RN=2*grad*diag.inc_uz.^2./(diag.inc_ur*cw^2);
    add([tag '_native_RN_incident_gradient'],max(abs(diag.wall_rn-expected_RN)),1e-16);
    pq_error=diag.p1_ref-diag.p1_in-diag.q1_in.*expected_RN;
    add([tag '_native_p_kick_error'],max(abs(pq_error)),1e-12);
    add([tag '_q_reflect_error'],max(diag.q_ref_error),1e-12);
    add([tag '_incident_unit_error'],max(abs(hypot(diag.inc_ur,diag.inc_uz)-1)),1e-7);
    coefficient_native=diag.amp_ref./diag.amp_in.*exp(1i*diag.phase_delta);
    omega=2*pi*cfg.frequency_hz;
    kappa=omega*diag.inc_uz/cw;
    physical_coefficient=reflection(kappa,cw,cfg);
    add([tag '_plane_wave_reflection_coefficient'],max(abs(conj(coefficient_native)-physical_coefficient)),2e-7);
    [~,ia]=min(abs(diag.alpha_deg));
    expected_time=integral(@(r)1./fold_c(r,cfg),0,max(cfg.receiver_range_m), ...
        'Waypoints',cfg.wall_range_m,'AbsTol',1e-12,'RelTol',1e-12);
    add([tag '_axis_time_error_s'],abs(diag.tau_receiver_real(ia)-expected_time),.05/omega);
    % Independent Snell shooting locates the stationary ray per requested point.
    p_ref=fields{2}.pressure; p_wall=fields{1}.pressure;
    expected_factor=zeros(size(p_ref)); cs=munk_profile_sound_speed(cfg.source_z_m,cfg.munk_parameters);
    for iz=1:size(p_ref,1)
        target=fields{1}.receiver_depth_m(iz)-cfg.transverse_origin_m;
        for ir=1:size(p_ref,2)
            R=cfg.receiver_range_m(ir);
            if target==0
                slowness=0;
            else
                bracket=sin(deg2rad(max(abs(cfg.angle_limits_deg))))/cs;
                fun=@(p)integral(@(r)p*fold_c(r,cfg)./sqrt(1-(p*fold_c(r,cfg)).^2),0,R, ...
                    'Waypoints',cfg.wall_range_m*(R>cfg.wall_range_m),'AbsTol',1e-10,'RelTol',1e-11)-target;
                slowness=fzero(fun,[-bracket bracket]);
            end
            expected_factor(iz,ir)=conj(reflection(omega*slowness,cw,cfg));
        end
    end
    expected=p_ref.*expected_factor;
    add([tag '_reflected_unfolded_complex_relative'],norm(p_wall(:)-expected(:))/norm(expected(:)),3e-4);
    add([tag '_reflected_boundary_samples_nonzero'],double(any(abs(p_wall(:,1))==0)),0);
    results{id}=struct('cfg',cfg,'metas',{metas},'fields',{fields},'diag',diag, ...
        'expected_RN',expected_RN,'coefficient_native',coefficient_native,'physical_coefficient',physical_coefficient, ...
        'expected_factor',expected_factor,'processes',{processes}, ...
        'note','Single-reflected component only; folded continuation is not independent FFP; native R spreading remains diagnostic.');
    save(fullfile(options.output_dir,[tag '_result.mat']),'cfg','meta','diag','p_wall','p_ref','expected_factor');
    fprintf('%s reflection complex relative=%.9g\n',tag,norm(p_wall(:)-expected(:))/norm(expected(:)));
end
status='COMPLETE';
catch exception
err=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'elapsed_s',toc(t),'result_dir',options.output_dir,'matlab_release',version, ...
    'overlay_exe_sha256',black_file_hash(options.overlay_exe), ...
    'scope','One physical flat wall; folded-Q internal consistency and plane-wave reflection. No complete total field or 3D spreading qualification.');
save(fullfile(options.output_dir,'result.mat'),'checks','receipt','results');
fid=fopen(fullfile(options.output_dir,'error.log'),'w'); fprintf(fid,'%s',err); fclose(fid);
fid=fopen(fullfile(options.output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
if ~isempty(err), fprintf('%s\n',err); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function c=fold_c(r,cfg)
sg=1-2*strcmp(cfg.direction,'up'); z=cfg.source_z_m+sg*r;
if strcmp(cfg.reflection_wall,'surface'), z=abs(z);
else, z=cfg.water_depth_m-abs(cfg.water_depth_m-z); end
c=munk_profile_sound_speed(z,cfg.munk_parameters);
end
function R=reflection(kappa,cw,cfg)
if strcmp(cfg.reflection_wall,'surface'), R=-ones(size(kappa)); return; end
omega=2*pi*cfg.frequency_hz;
alpha=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
kb=omega/cfg.bottom_cp_mps+1i*alpha;
qw=sqrt((omega/cw)^2-kappa.^2); qb=sqrt(kb^2-kappa.^2);
qb(imag(qb)<0)=-qb(imag(qb)<0);
R=(qw-qb/cfg.bottom_rho_gcm3)./(qw+qb/cfg.bottom_rho_gcm3);
end