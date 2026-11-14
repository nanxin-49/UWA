function receipt = validate_bellhop_rotated_multifield(options)
%VALIDATE_BELLHOP_ROTATED_MULTIFIELD Sparse explicit-order/phase/carry gate.
% Orders 10/12 are coarse-Q cutoff/count controls, not final Munk accuracy.
arguments
    options.overlay_exe (1,:) char
    options.output_dir (1,:) char
    options.case_indices (1,:) double {mustBeInteger,mustBePositive} = 1:8
    options.long_step_m (1,1) double {mustBePositive} = 2
end
assert(~isfolder(options.output_dir)); mkdir(options.output_dir); t=tic;
diary(fullfile(options.output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
directions={'up','down','up','down','up','down','up','down'};
orders=[1 1 2 2 3 3 10 12];
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
runs=cell(1,numel(orders)); status='FAILED'; err='';
try
assert(all(options.case_indices<=numel(orders)));
for j=options.case_indices
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM gate.');
    assert(toc(t)<60,'Multiwall sparse pilot cumulative time gate.');
    cfg=make_bellhop_rotated_multifield_config(directions{j},orders(j),[0 1000 3000 4900 5000],[-1 0 1]);
    if orders(j)>=10, cfg.ssp_spacing_m=20; cfg.step_m=options.long_step_m; cfg.beam_count=101; cfg.angle_limits_deg=[-.02 .02]; end
    tag=sprintf('%s_order%d',directions{j},orders(j)); fields=cell(1,2); metas=cell(1,2); processes=cell(1,2);
    for ib=1:2
        folder=fullfile(options.output_dir,sprintf('%s_mode%d',tag,ib)); mkdir(folder); stem=fullfile(folder,'case');
        metas{ib}=write_bellhop_rotated_multifield_q(stem,cfg);
        if ib==2
            fid=fopen([stem '.iwcov'],'w'); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
        end
        processes{ib}=run_bellhop_rotated_local_case(options.overlay_exe,stem,min(20,60-toc(t)));
        fields{ib}=read_bellhop_rotated_field(stem);
    end
    folder=fullfile(options.output_dir,[tag '_mode1']); stem=fullfile(folder,'case');
    mw=readmatrix([stem '.mwdiag'],'FileType','text','CommentStyle','#'); mw=mw(~all(isnan(mw),2),:);
    d=readmatrix([stem '.r25diag'],'FileType','text','CommentStyle','#'); d=d(~all(isnan(d),2),:);
    assert(height(mw)==cfg.beam_count && size(mw,2)==8 && size(d,2)==10);
    add([tag '_requested_count'],max(abs(mw(:,2)-cfg.reflection_order)),0);
    add([tag '_walls_completed'],max(abs(mw(:,3)-cfg.reflection_order)),0);
    add([tag '_physical_surface_count'],max(abs(mw(:,4)-cfg.expected_surface_bounces)),0);
    add([tag '_physical_bottom_count'],max(abs(mw(:,5)-cfg.expected_bottom_bounces)),0);
    add([tag '_native_counter_agreement'],max(abs(d(:,9:10)-mw(:,4:5)),[],'all'),0);
    add([tag '_trajectory_capacity_fraction'],max(mw(:,6))/600000,.95);
    add([tag '_receiver_range_reached'],double(any(mw(:,7)<max(cfg.receiver_range_m))),0);
    parity=(-1)^cfg.reflection_order;
    p=sind(d(:,1))./d(:,7);
    carry=(d(:,3)-cfg.transverse_origin_m)-parity*p.*d(:,6);
    add([tag '_outplane_carry_relative'],max(abs(carry))/max(1,max(abs(d(:,3)-cfg.transverse_origin_m))),1e-8);
    expected_factor=bellhop_rotated_multifield_reference_factor(cfg,metas{1}, ...
        cfg.receiver_range_m,fields{1}.receiver_depth_m);
    expected=fields{2}.pressure.*expected_factor;
    relative=norm(fields{1}.pressure(:)-expected(:))/norm(expected(:));
    add([tag '_reflected_unfolded_field_relative'],relative,3e-4);
    add([tag '_boundary_fields_nonzero'],double(any(abs(fields{1}.pressure(:,[1 end]))==0,'all')),0);
    if cfg.reflection_order>=10
        add([tag '_below_legacy_amp_cutoff'],double(any(mw(:,8)>=.005)),0);
    end
    pressure_ffp=-conj(fields{1}.pressure);
    physical_pressure=complex(zeros(numel(cfg.nominal_receiver_depths_m),numel(cfg.nominal_receiver_x_m)));
    physical_pressure(cfg.receiver_depth_order,cfg.receiver_x_order)=pressure_ffp.';
    runs{j}=struct('cfg',cfg,'metas',{metas},'fields',{fields},'processes',{processes}, ...
        'wall_diagnostics',mw,'ray_diagnostics',d,'expected_factor',expected_factor, ...
        'physical_pressure',physical_pressure,'relative',relative, ...
        'scope','Sparse component only; coarse order10/12 cases do not qualify Munk phase accuracy');
    save(fullfile(options.output_dir,[tag '_result.mat']),'cfg','physical_pressure','mw','d','expected_factor','relative');
    fprintf('%s field_relative=%.9g steps_max=%d amp_max=%.9g\n',tag,relative,max(mw(:,6)),max(mw(:,8)));
end
status='COMPLETE';
catch exception
err=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'elapsed_s',toc(t),'result_dir',options.output_dir,'matlab_release',version, ...
    'case_indices',options.case_indices,'long_step_m',options.long_step_m, ...
    'scope','Sparse explicit reflection count/carry/coefficient controls; no complete field', ...
    'overlay_exe_sha256',black_file_hash(options.overlay_exe));
save(fullfile(options.output_dir,'result.mat'),'checks','receipt','runs');
fid=fopen(fullfile(options.output_dir,'error.log'),'w'); fprintf(fid,'%s',err); fclose(fid);
fid=fopen(fullfile(options.output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt)); if ~isempty(err), fprintf('%s\n',err); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
