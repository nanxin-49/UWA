function receipt = validate_bellhop_rotated_source_charts(options)
%VALIDATE_BELLHOP_ROTATED_SOURCE_CHARTS Direct source plane and chart overlaps.
arguments
    options.overlay_exe (1,:) char
    options.output_dir (1,:) char
    options.reuse_dir (1,:) char = ''
end
assert(~isfolder(options.output_dir)); mkdir(options.output_dir); t=tic;
diary(fullfile(options.output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
runs=cell(1,6); err=''; status='FAILED';
R=[.05 .1 1 5]; Q=[-5 -1 -.1 -.05 0 .05 .1 1 5];
tags={'uniform_right','munk_right_native','munk_right_point','munk_left_point','munk_up_point','munk_down_point'};
try
for j=1:6
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM gate.');
    assert(toc(t)<40,'Source chart cumulative time gate.');
    cfg=make_bellhop_rotated_munk_config('up'); cfg.run_type='C';
    cfg.uniform_control=j==1; cfg.point3d=j~=2;
    cfg.beam_count=2001; cfg.angle_limits_deg=[-46 46]; cfg.step_m=.02;
    cfg.range_limit_m=5.1; cfg.receiver_range_m=R; cfg.receiver_q_m=Q;
    folder=fullfile(options.output_dir,tags{j}); mkdir(folder); stem=fullfile(folder,'case');
    if j<=4
        side=1; if j==4, side=-1; end
        meta=write_bellhop_rotated_horizontal_q(stem,cfg,side);
    else
        cfg.direction='up'; if j==6, cfg.direction='down'; end
        cfg.receiver_x_m=Q; if j==6, cfg.receiver_x_m=-Q; end
        meta=write_bellhop_rotated_munk_q(stem,cfg);
        fid=fopen([stem '.iwcov'],'w'); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
        write_bellhop_rotated_field_options(stem,R,meta.q_receiver_m,point3d=true,fast_q=true);
    end
    reuse=false;
    if ~isempty(options.reuse_dir)
        prior_cp=jsondecode(fileread(fullfile(options.reuse_dir,'checkpoint.json')));
        assert(strcmp(prior_cp.overlay_exe_sha256,black_file_hash(options.overlay_exe)));
        prior=fullfile(options.reuse_dir,tags{j});
        if isfile(fullfile(prior,'case.shd')) && isfile(fullfile(prior,'state.txt'))
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
    end
    if ~reuse
        process=run_bellhop_rotated_local_case(options.overlay_exe,stem,min(20,40-toc(t)));
    end
    field=read_bellhop_rotated_field(stem);
    d=readmatrix([stem '.r25diag'],'FileType','text','CommentStyle','#');
    d=d(~all(isnan(d),2),:);
    add([tags{j} '_artificial_bounce_count'],max(sum(d(:,9:10),2)),0);
    if j<=4
        identity=cfg.munk_parameters.c0_mps; %#ok<NASGU>
        expected_I=d(:,7).*d(:,2)./cosd(d(:,1));
        add([tags{j} '_natural_chart_outplane_relative'],max(abs(d(:,6)-expected_I)./expected_I),1e-9);
    end
    coverage=abs(Q(:))<=R;
    raw_pressure=field.pressure; covered_pressure=raw_pressure; covered_pressure(~coverage)=NaN;
    runs{j}=struct('cfg',cfg,'meta',meta,'field',field,'process',process, ...
        'coverage',coverage,'covered_pressure',covered_pressure,'ray_diagnostics',d, ...
        'note','Raw pressure retained. NaN marks points outside this chart sector, not missing-axis interpolation.');
    save(fullfile(folder,'result.mat'),'cfg','meta','field','process','coverage','covered_pressure','raw_pressure','d');
end
distance=hypot(Q(:),R);
exact=exp(1i*2*pi*10000*distance/1500)./distance;
p=-conj(runs{1}.field.pressure); use=runs{1}.coverage;
add('uniform_source_chart_Green_complex',norm(p(use)-exact(use))/norm(exact(use)),1e-3);
add('natural_horizontal_point_vs_native', ...
    norm(runs{2}.field.pressure(use)-runs{3}.field.pressure(use))/norm(runs{2}.field.pressure(use)),1e-4);
left=flipud(runs{4}.field.pressure); right=runs{3}.field.pressure;
add('left_right_proper_rotation_symmetry',norm(left(use)-right(use))/norm(right(use)),1e-4);
for ir=1:numel(R)
    ip=find(abs(Q-R(ir))<1e-12); im=find(abs(Q+R(ir))<1e-12);
    pairs=[runs{5}.field.pressure(ip,ir),runs{3}.field.pressure(im,ir); ...
        runs{5}.field.pressure(im,ir),runs{4}.field.pressure(ip,ir); ...
        runs{6}.field.pressure(im,ir),runs{3}.field.pressure(ip,ir); ...
        runs{6}.field.pressure(ip,ir),runs{4}.field.pressure(im,ir)];
    add(sprintf('overlap_at_%.2gm_relative',R(ir)),max(abs(pairs(:,1)-pairs(:,2))./abs(pairs(:,2))),1e-3);
end
row=find(Q==0);
add('source_depth_plane_all_non_source_points_finite',double(any(~isfinite([ ...
    runs{3}.field.pressure(row,:),runs{4}.field.pressure(row,:)]))),0);
status='COMPLETE';
catch exception
err=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'elapsed_s',toc(t),'result_dir',options.output_dir,'matlab_release',version, ...
    'scope','Small source-plane and chart-overlap qualification; no reflected/complete field', ...
    'overlay_exe_sha256',black_file_hash(options.overlay_exe));
save(fullfile(options.output_dir,'result.mat'),'receipt','checks','runs');
fid=fopen(fullfile(options.output_dir,'error.log'),'w'); fprintf(fid,'%s',err); fclose(fid);
fid=fopen(fullfile(options.output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt)); if ~isempty(err), fprintf('%s\n',err); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
