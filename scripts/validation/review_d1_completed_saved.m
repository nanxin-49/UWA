function receipt=review_d1_completed_saved(origin,tx100_prior,out)
%REVIEW_D1_COMPLETED_SAVED Final independent saved-data and local FIG closeout.
arguments
    origin (1,:) char
    tx100_prior (1,:) char
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out); diary(fullfile(out,'diary.log')); cleanup=onCleanup(@()diary('off')); t0=tic;
previous_threads=maxNumCompThreads; maxNumCompThreads(1); restore_threads=onCleanup(@()maxNumCompThreads(previous_threads));
rows=jsondecode(fileread(fullfile(origin,'retrieval','manifest.json'))); pm_sensitivity_io('verify',origin,rows);
cp=jsondecode(fileread(fullfile(origin,'output','checkpoint.json')));
assert(strcmp(cp.status,'D1A_NUMERICS_COMPLETE_LOCAL_FIGURES_REQUIRED') && cp.new_BIE_solves==5 && cp.imported_BIE_solves==5);
assert(numel(cp.completed_steps)==12 && isscalar(cp.failed_or_unrun_checks));
parent=jsondecode(fileread(fullfile(tx100_prior,'checkpoint.json'))); pm_sensitivity_io('verify',tx100_prior,parent.artifacts_sha256);
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'}; recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96];
all_metrics=table; physical_change=table; figure_curves=0; scales=zeros(2,1); Q=cell(2,1); A=Q;
for index=1:2
    depths=[100 1490]; depth=depths(index); src=fullfile(origin,'output',sprintf('depth_%04dm',depth));
    work=fullfile(out,sprintf('depth_%04dm',depth)); mkdir(work);
    for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','policy.json'}
        copyfile(fullfile(src,name{1}),fullfile(work,name{1}));
    end
    copyfile(fullfile(src,'input'),fullfile(work,'input')); copyfile(fullfile(src,'jobs'),fullfile(work,'jobs'));
    q=load(fullfile(work,'common.mat')); a=load(fullfile(work,'comparison.mat')); Q{index}=q; A{index}=a;
    assert(q.cfg.source.z_source_m==depth && q.cfg.source.total_width_m==.1 && q.cfg.seed==260001 && all(q.z==6));
    assert(numel(q.x)==411 && max(abs(q.x-linspace(-36,36,411)'))<1e-12 && a.reference_qualified);
    if depth==100
        % Prior local audit is reusable only when every mathematical input is exact.
        for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv'}
            assert(strcmp(pm_sensitivity_io('sha',fullfile(src,name{1})),pm_sensitivity_io('sha',fullfile(tx100_prior,'audited_work',name{1}))));
        end
        for j=1:5
            for name={'data.mat','points.csv'}
                assert(strcmp(pm_sensitivity_io('sha',fullfile(src,'jobs',ids{j},name{1})),pm_sensitivity_io('sha',fullfile(tx100_prior,'audited_work','jobs',ids{j},name{1}))));
            end
        end
        mode='REUSED_IDENTICAL_INPUTS_LOCAL_NUMERICAL_AUDIT';
    else
        for j=1:5
            stage=fullfile(work,'jobs',ids{j}); seal=jsondecode(fileread(fullfile(stage,'checkpoint.json')));
            assert(strcmp(seal.status,'AUDITED_PASS')); pm_sensitivity_io('verify',stage,seal.artifacts_sha256);
            audit_r1_bie_manual_stage(stage,fullfile(work,'common.mat'),recipes(j,:));
        end
        audit_d1_models(work); diary(fullfile(out,'diary.log')); mode='NEW_LOCAL_SAVED_DENSITY_AND_MODEL_AUDIT';
    end
    audit_d1_final(work); render_d1_depth_results(work); diary(fullfile(out,'diary.log'));
    fa=jsondecode(fileread(fullfile(work,'figures','figure_audit.json'))); assert(fa.failed_checks==0 && fa.differing_pixel_channels==0);
    figure_curves=figure_curves+fa.reopened_curves;
    M=readtable(fullfile(work,'model_metrics.csv')); scales(index)=norm(a.U(:,2))/sqrt(411);
    assert(abs(scales(index)-a.scale)<1e-12 && isequal(a.mask,abs(a.U(:,2))>=.01*scales(index)));
    rough_flat=norm(a.U(:,2)-q.flat)/norm(q.flat); assert(abs(rough_flat-a.rough_flat)<1e-12);
    physical_change=[physical_change;table(depth,rough_flat,scales(index),'VariableNames',{'tx_depth_m','BIE_rough_flat_relative_L2','BIE_RMS_Pa'})]; %#ok<AGROW>
    for method=[1 3]
        error_f=a.U(:,method)-a.U(:,2); amp=abs(a.U(:,method))-abs(a.U(:,2)); ph=angle(a.U(:,method).*conj(a.U(:,2)));
        E=norm(error_f)/norm(a.U(:,2)); phase_component=norm(2*sqrt(abs(a.U(:,method)).*abs(a.U(:,2))).*sin(ph/2))/norm(a.U(:,2));
        assert(abs(E-M.complex_relative_L2(method))<1e-12 && abs(norm(amp)/norm(a.U(:,2))-M.amplitude_relative_L2(method))<1e-12);
        assert(abs(sqrt(mean(ph(a.mask).^2))-M.phase_rms_rad(method))<1e-12);
        assert(abs(E^2-M.amplitude_relative_L2(method)^2-phase_component^2)<1e-11);
        row=M(method,:); row.tx_depth_m=depth; row.BIE_RMS_Pa=scales(index); all_metrics=[all_metrics;row]; %#ok<AGROW>
    end
    pm_sensitivity_io('write',fullfile(work,'local_review.json'),struct('status','PASS','failed_checks',0,'mode',mode,'new_BIE_solves',0,'figure_audit',fa));
end
assert(isequal(Q{1}.profile,Q{2}.profile) && isequal(Q{1}.x,Q{2}.x) && isequal(Q{1}.z,Q{2}.z));
N=readtable(fullfile(origin,'output','depth_comparison.csv'));
for j=1:height(all_metrics)
    ii=string(N.method)==string(all_metrics.method(j)) & N.tx_depth_m==all_metrics.tx_depth_m(j);
    assert(nnz(ii)==1 && abs(N.complex_relative_L2(ii)-all_metrics.complex_relative_L2(j))<1e-12);
end
D=readtable(fullfile(origin,'output','depth_change.csv'));
for j=1:2
    v=all_metrics.complex_relative_L2(j+2)-all_metrics.complex_relative_L2(j);
    limit=all_metrics.five_times_empirical_resolution(j+2)+all_metrics.five_times_empirical_resolution(j);
    assert(abs(D.E1490_minus_E100(j)-v)<1e-12 && abs(D.conservative_change_resolution(j)-limit)<1e-12 && logical(D.change_resolved(j))==(abs(v)>limit));
end
write_d1_csv(all_metrics,fullfile(out,'depth_summary_review.csv')); write_d1_csv(D,fullfile(out,'depth_change_review.csv'));
write_d1_csv(physical_change,fullfile(out,'physical_change_review.csv'));
save(fullfile(out,'summary.mat'),'all_metrics','D','scales','physical_change','-v7');
pm_sensitivity_io('verify',origin,rows); pm_sensitivity_io('verify',tx100_prior,parent.artifacts_sha256);
receipt=struct('status','D1A_LOCAL_NUMERICAL_AND_FIGURE_AUDIT_PASS','failed_checks',0,'elapsed_s',toc(t0), ...
    'origin_SHA_files',numel(rows),'new_BIE_solves',0,'independent_reopened_curves',figure_curves,'different_pixel_channels',0, ...
    'matlab_version',version,'remaining_required_checks',{{}},'scope','one paired original PM realization; fixed 2D 4kHz Rx6 line');
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt);
files=dir(fullfile(out,'**','*')); paths={};
for j=1:numel(files)
    if ~files(j).isdir && ~strcmp(files(j).name,'diary.log') && ~strcmp(files(j).name,'checkpoint.json')
        p=fullfile(files(j).folder,files(j).name); paths{end+1}=strrep(p(numel(out)+2:end),'\','/'); %#ok<AGROW>
    end
end
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
pm_sensitivity_io('write',fullfile(out,'checkpoint.json'),struct('status','AUDITED_PASS','artifacts_sha256',pm_sensitivity_io('manifest',out,paths), ...
    'code_versions_sha256',pm_sensitivity_io('manifest',repo,{'scripts/validation/review_d1_completed_saved.m','scripts/validation/audit_d1_models.m','scripts/reporting/render_d1_depth_results.m'}), ...
    'failed_or_unrun_checks',{{}},'new_BIE_solves',0,'next_allowed_step','Write final D1 report; D1b remains feasibility only'));
fprintf('%s\n',jsonencode(receipt));
end
