function receipt = validate_bellhop_rotated_munk_rays(options)
%VALIDATE_BELLHOP_ROTATED_MUNK_RAYS Independent Snell/integral path audit.
% Uses three rays only and compares each retained node directly, no field fill.
arguments
    options.bellhop_exe (1,:) char
    options.output_dir (1,:) char
end
assert(~isfolder(options.output_dir)); mkdir(options.output_dir);
diary(fullfile(options.output_dir,'diary.log')); diary_cleanup=onCleanup(@()diary('off'));
t=tic; checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
runs=cell(1,2); directions={'up','down'}; err=''; status='FAILED';
try
for id=1:2
    cfg=make_bellhop_rotated_munk_config(directions{id});
    cfg.run_type='R'; cfg.beam_count=3; cfg.step_m=.5; cfg.ssp_spacing_m=2.5;
    if id==1, cfg.angle_limits_deg=[-.05 .05]; else, cfg.angle_limits_deg=[-5 5]; end
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM preflight gate.');
    folder=fullfile(options.output_dir,directions{id}); mkdir(folder);
    case_root=fullfile(folder,'case'); meta=write_bellhop_rotated_munk_q(case_root,cfg);
    proc=run_bellhop_rotated_local_case(options.bellhop_exe,case_root,20);
    rays=read_rays([case_root '.ray']);
    add([directions{id} '_three_rays'],abs(numel(rays)-3),0);
    c_source=munk_profile_sound_speed(cfg.source_z_m,cfg.munk_parameters);
    max_path=0; max_invariant=0; max_boundary=0; max_reversal=0; max_map=0;
    for j=1:numel(rays)
        a=deg2rad(rays(j).angle_deg);
        p=sin(a)/c_source; % conserved transverse slowness in c(r)
        r=rays(j).rq(:,1); q=rays(j).rq(:,2)-cfg.transverse_origin_m;
        keep=r>=0 & r<=cfg.range_limit_m;
        r=r(keep); q=q(keep);
        sg=1-2*strcmp(cfg.direction,'up');
        z=cfg.source_z_m+sg*r; physical_x=-sg*q;
        % Integrate the continuous Munk profile independently of the Q table.
        x_exact=arrayfun(@(R)integral(@(rr)integrand(rr,p,sg,cfg),0,R, ...
            'RelTol',1e-11,'AbsTol',1e-11),r);
        max_path=max(max_path,max(abs(q-x_exact)));
        % Output coordinates are ASCII single precision. Chord slopes are
        % therefore tested only over segments with >=0.1m range separation.
        dr=diff(r); dq=diff(q); good=dr>=.1;
        mid=.5*(r(1:end-1)+r(2:end));
        cm=munk_profile_sound_speed(cfg.source_z_m+sg*mid,cfg.munk_parameters);
        observed=dq./hypot(dr,dq)./cm;
        max_invariant=max(max_invariant,max(abs(observed(good)-p)));
        max_boundary=max(max_boundary,max(abs(q))/cfg.artificial_halfwidth_m);
        max_reversal=max(max_reversal,double(any(dr<=0)));
        % Inverse transform and forward transform without interpolating an axis.
        r2=sg*(z-cfg.source_z_m); q2=-sg*physical_x;
        max_map=max(max_map,max(hypot(r-r2,q-q2)));
        assert(rays(j).top==0 && rays(j).bottom==0,'Artificial boundary bounce.');
    end
    add([directions{id} '_path_error_m'],max_path,2e-3);
    add([directions{id} '_transverse_slowness_error_s_per_m'],max_invariant,2e-6);
    add([directions{id} '_artificial_boundary_fraction'],max_boundary,.1);
    add([directions{id} '_range_reversal_count'],max_reversal,0);
    add([directions{id} '_inverse_map_error_m'],max_map,1e-11);
    runs{id}=struct('cfg',cfg,'meta',meta,'rays',rays,'process',proc);
end
status='COMPLETE';
catch exception
    err=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
fid=fopen(fullfile(options.output_dir,'error.log'),'w'); fprintf(fid,'%s',err); fclose(fid);
receipt=struct('status',status,'elapsed_s',toc(t),'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'scope','Three direct rays in each rotated chart; not field or reflection qualification', ...
    'result_dir',options.output_dir,'matlab_release',version, ...
    'audit_note','Continuous Snell quadrature independent of Q tabulation; ASCII coordinate precision retained');
save(fullfile(options.output_dir,'result.mat'),'receipt','checks','runs');
fid=fopen(fullfile(options.output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function f=integrand(r,p,sg,cfg)
c=munk_profile_sound_speed(cfg.source_z_m+sg*r,cfg.munk_parameters);
f=p*c./sqrt(1-(p*c).^2);
end
function rays=read_rays(file)
fid=fopen(file,'r'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fgetl(fid); fscanf(fid,'%f',1); fscanf(fid,'%d',3); beams=fscanf(fid,'%d',2);
fscanf(fid,'%f',2); fgetl(fid); % consume remaining header newline
kind=strtrim(fgetl(fid)); assert(contains(kind,'rz'),'Expected rz ray output.');
rays=struct('angle_deg',{},'top',{},'bottom',{},'rq',{});
for j=1:beams(1)
    angle=fscanf(fid,'%f',1); info=fscanf(fid,'%d',3);
    assert(numel(info)==3,'Incomplete ray header.');
    rq=fscanf(fid,'%f',[2 info(1)]).';
    rays(j)=struct('angle_deg',angle,'top',info(2),'bottom',info(3),'rq',rq);
end
end