function receipt=audit_munk_ffp_full(run_id,audit_id)
%AUDIT_MUNK_FFP_FULL Re-open saved fields and independently audit requirements.
arguments
    run_id (1,:) char = '20261112_full_01'
    audit_id (1,:) char = 'audit_01'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'scripts','validation','support'),'-begin');
assert(isempty(regexp([run_id audit_id],'[^A-Za-z0-9_-]','once')));
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id); a=fullfile(d,audit_id);
assert(~isfolder(a),'Preserve every audit attempt'); mkdir(a); diary(fullfile(a,'diary.log')); t=tic;
receipt=struct('stage','full_saved_audit','status','RUNNING','elapsed_s',0,'exit_status',NaN, ...
    'qualification','UNVERIFIED','failed_checks',0,'failed_or_unrun_checks',{{}});
try
    has_refined=isfile(fullfile(d,'J07.mat')); final_job=6+has_refined;
    files=arrayfun(@(j)fullfile(d,sprintf('J%02d.mat',j)),1:final_job,'UniformOutput',false);
    assert(all(cellfun(@isfile,files)),'All six planned numerical fields are required');
    frozen=load(fullfile(d,'frozen_plan.mat'),'plan');
    record=struct([]); comparisons=[1 2;2 3;3 4;4 5;5 6];
    names={'early_layer','cutoff','quadrature','contour','final_layer'};
    if has_refined
        refinement=load(fullfile(d,'refinement_plan.mat'),'cfg');
        comparisons(end+1,:)=[6 7]; names{3}='early_quadrature'; names{6}='final_quadrature';
    end
    for pair=1:size(comparisons,1)
        p=local_read(files{comparisons(pair,1)}); q=local_read(files{comparisons(pair,2)});
        assert(isequal(p.z,q.z) && isequal(p.x,q.x));
        expected_a=frozen.plan(comparisons(pair,1)).config;
        if comparisons(pair,2)==7, expected_b=refinement.cfg; else, expected_b=frozen.plan(comparisons(pair,2)).config; end
        assert(isequaln(p.cfg,expected_a) && isequaln(q.cfg,expected_b),'Saved config differs from frozen plan');
        masks={p.z<4900 & abs(p.x)>=0,p.z>=4900 & abs(p.x)>=0, ...
            true(size(p.z)) & p.x==0,true(size(p.z)) & abs(p.x)>=.05};
        region_names={'upper_column','near_bottom','axis','off_axis'};
        for j=1:4
            mask=masks{j} & isfinite(p.g) & isfinite(q.g);
            pp=p.g(mask); qq=q.g(mask); peak=max(abs(qq));
            rel=sqrt(sum(abs(pp-qq).^2)/sum(abs(qq).^2));
            good=abs(pp)>peak*1e-6 & abs(qq)>peak*1e-6;
            phase=abs(atan2(imag(pp(good).*conj(qq(good))),real(pp(good).*conj(qq(good)))));
            db=abs(20*log10(abs(pp(good))./abs(qq(good))));
            pass=rel<=1e-3 && local_q(phase,.95)<=.01 && local_q(db,.95)<=.05;
            record=[record struct('comparison',string(names{pair}),'a_job',comparisons(pair,1), ...
                'b_job',comparisons(pair,2),'region',string(region_names{j}),'complex_rel_l2',rel, ...
                'max_abs_green_difference',max(abs(pp-qq)),'phase_p95_rad',local_q(phase,.95), ...
                'tl_p95_db',local_q(db,.95),'valid_samples',nnz(mask),'phase_tl_support_fraction',mean(good), ...
                'required_final_gate',pair~=1 && (~has_refined || pair~=3),'pass',pass)]; %#ok<AGROW>
        end
        clear p q;
    end
    tbl=struct2table(record); writetable(tbl,fullfile(a,'convergence.csv'));
    s=load(files{final_job},'output'); o=s.output; g=o.green; mask=~o.source_singularity_mask;
    assert(isequal(size(g),[10981 201]));
    assert(nnz(o.source_singularity_mask)==1 && all(isfinite(g(mask))) && isnan(g(o.source_singularity_mask)));
    scale_error=max(abs(o.pressure_1m(mask)-4*pi*g(mask)));
    axis=o.range_m==0; valid_axis=~o.source_singularity_mask(:,axis);
    axis_finite=all(isfinite(g(valid_axis,axis)));
    mirror=fliplr(g); parity=sqrt(sum(abs(g(mask)-mirror(mask)).^2)/sum(abs(g(mask)).^2));
    surface=max(abs(g(o.depth_m==0,:)),[],'all')/max(abs(g(mask)));
    localz=o.depth_m(o.depth_m>=4900); maxdz=max(diff(localz));
    lambda_min=min(munk_profile_sound_speed(localz,o.cfg.munk_parameters))/o.cfg.frequency_hz;
    tau_ref=munk_vertical_travel_time_reference([0 1000 3000 4900 5000],4950,o.cfg.munk_parameters);
    tau_model=zeros(size(tau_ref)); testz=[0 1000 3000 4900 5000];
    for j=1:numel(testz)
        lo=min(testz(j),4950); hi=max(testz(j),4950);
        lengths=max(0,min(o.cfg.layer_edges_m(2:end),hi)-max(o.cfg.layer_edges_m(1:end-1),lo));
        tau_model(j)=sum(lengths./o.cfg.layer_c_mps);
    end
    tau_error=max(abs(tau_model-tau_ref));
    w=2*pi*o.cfg.frequency_hz;
    cumulative=[0;cumsum(diff(o.cfg.layer_edges_m)./o.cfg.layer_c_mps)];
    tau_axis=abs(interp1(o.cfg.layer_edges_m,cumulative,o.depth_m)-interp1(o.cfg.layer_edges_m,cumulative,4950));
    paxis=o.pressure_1m(:,axis); phase=angle(paxis); phase(~isfinite(paxis))=NaN;
    axistable=table(o.depth_m,real(paxis),imag(paxis),abs(paxis),phase, ...
        angle(paxis.*exp(-1i*w*tau_axis)),o.source_singularity_mask(:,axis), ...
        'VariableNames',{'z_m','real_pressure_1m','imag_pressure_1m','magnitude','wrapped_phase_rad','onepass_referenced_phase_rad','true_source_singularity'});
    writetable(axistable,fullfile(a,'axis_complex.csv'));
    r=readtable(fullfile(a,'axis_complex.csv')); csv_error=max(abs(complex(r.real_pressure_1m(valid_axis),r.imag_pressure_1m(valid_axis))-paxis(valid_axis)));
    pre=fullfile(root,'results','validation','munk_ffp_5000m','20261112_preflight_01');
    sf=load(fullfile(pre,'uniform_free.mat'),'uniform_free'); sr=load(fullfile(pre,'uniform_release.mat'),'uniform_release');
    u=sf.uniform_free; [zz,xx]=ndgrid(u.depth_m,u.range_m); pts=[xx(:),zeros(numel(xx),1),zz(:)];
    exact=reshape(deepsea3d_point_field(pts,[0 0 12],2*pi*1e4/1500),size(xx));
    exact_mirror=complex(zeros(size(exact)));
    for j=1:numel(exact)
        pc=struct('frequency_hz',1e4,'c0_mps',1500,'tx_m',[0 0 12],'rx_m',pts(j,:));
        temp=deepsea3d_flat(pc,0); exact_mirror(j)=temp.total_native;
    end
    free_rel=norm(u.green-exact,'fro')/norm(exact,'fro');
    mirror_rel=norm(sr.uniform_release.green-exact_mirror,'fro')/norm(exact_mirror,'fro');
    free_tl=max(abs(20*log10(abs(u.green)./abs(exact))),[],'all');
    mirror_tl=max(abs(20*log10(abs(sr.uniform_release.green)./abs(exact_mirror))),[],'all');
    free_phase=max(abs(angle(u.green.*conj(exact))),[],'all');
    mirror_phase=max(abs(angle(sr.uniform_release.green.*conj(exact_mirror))),[],'all');
    long_reference=jsondecode(fileread(fullfile(d,'long_uniform_receipt.json')));
    assert(long_reference.exit_status==0 && long_reference.metrics.source_depth_m==4950 && ~long_reference.metrics.analytic_subtraction_used);
    gates=struct('uniform_actual_5000m_scale',long_reference.all_depths_pass, ...
        'complete_finite_axis',axis_finite,'all_non_source_receivers_finite',all(isfinite(g(mask))), ...
        'surface_boundary',surface<=1e-10,'left_right_parity',parity<=1e-12, ...
        'normalization_identity',scale_error==0,'high_resolution_dz',maxdz<=lambda_min/10, ...
        'vertical_onepass_time',w*tau_error<=.05,'uniform_free',free_rel<=1e-3 && free_tl<=.05 && free_phase<=.01, ...
        'uniform_mirror',mirror_rel<=1e-3 && mirror_tl<=.05 && mirror_phase<=.01, ...
        'final_cutoff_quadrature_contour_layer',all(tbl.pass(tbl.required_final_gate)), ...
        'axis_csv_roundtrip',csv_error<=1e-10);
    required=cell2mat(struct2cell(gates)); failed=fieldnames(gates); failed=failed(~required(:));
    metrics=struct('surface_rel',surface,'complex_parity_rel',parity,'max_patch_dz_m',maxdz,'lambda_min_m',lambda_min, ...
        'normalization_max_abs',scale_error,'axis_csv_max_abs',csv_error,'max_onepass_delay_error_s',tau_error, ...
        'uniform_free_rel',free_rel,'uniform_mirror_rel',mirror_rel,'uniform_free_max_tl_db',free_tl, ...
        'uniform_mirror_max_tl_db',mirror_tl,'uniform_free_max_phase_rad',free_phase,'uniform_mirror_max_phase_rad',mirror_phase);
    save(fullfile(a,'audit_values.mat'),'tbl','gates','metrics','tau_ref','tau_model','axistable','-v7');
    receipt.failed_checks=numel(failed); receipt.failed_or_unrun_checks=failed; receipt.gates=gates; receipt.metrics=metrics;
    if isempty(failed), receipt.status='REQUIRED_NUMERICAL_GATES_PASS'; receipt.qualification='CONDITIONAL_MODEL_NUMERICALLY_QUALIFIED'; receipt.exit_status=0;
    else, receipt.status='NUMERICAL_GATES_FAILED'; receipt.qualification='UNQUALIFIED'; receipt.exit_status=1; end
catch err
    receipt.status='AUDIT_FAILED'; receipt.exit_status=1; receipt.failed_checks=receipt.failed_checks+1;
    receipt.failed_or_unrun_checks={err.identifier,err.message};
    fid=fopen(fullfile(a,'error.log'),'w'); fprintf(fid,'%s',getReport(err,'extended','hyperlinks','off')); fclose(fid);
end
if exist('final_job','var'), receipt.final_job=final_job; else, receipt.final_job=NaN; end; receipt.elapsed_s=toc(t); receipt.audit_dir=a; local_json(fullfile(a,'receipt.json'),receipt); diary off;
fprintf('%s\n',jsonencode(receipt));
end
function p=local_read(path)
s=load(path,'output'); p=struct('g',s.output.green,'z',s.output.depth_m,'x',s.output.range_m,'cfg',s.output.cfg);
end
function value=local_q(x,q)
x=sort(x(:)); if isempty(x), value=NaN; else, value=x(ceil(q*numel(x))); end
end
function local_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
