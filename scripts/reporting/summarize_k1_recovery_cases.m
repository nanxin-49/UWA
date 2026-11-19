function receipt=summarize_k1_recovery_cases(master,out,mode,qualification_file)
%SUMMARIZE_K1_PM_CASES Audit-bound seven-condition report from saved fields.
if nargin<3, mode='real'; end
if nargin<4, qualification_file=''; end
assert(ismember(mode,{'real','test-only','local-closeout'})); test_only=strcmp(mode,'test-only');
local_closeout=strcmp(mode,'local-closeout');
assert(~isfolder(out)); mkdir(out); timer=tic;
manifest=k1_pm_io('read',fullfile(master,'manifest.json')); pm_sensitivity_io('verify',master,manifest.files);
cfg=k1_pm_io('read',fullfile(master,'config.json')); proof=k1_pm_io('read',fullfile(master,'anchor','anchor_proof.json'));
if test_only, assert(strcmp(cfg.kind,'TEST_ONLY_K1_SUMMARY'));
else, assert(strcmp(cfg.kind,'K1_PM_RECOVERY_FIVE_REMAINING')); end
assert(strcmp(proof.status,'STRICT_R1_K0_ANCHOR_MATCH_PASS'));
if local_closeout
    qualifications=k1_pm_io('read',qualification_file);
    assert(strcmp(qualifications.status,'K1_SIX_CASES_INDEPENDENT_CLOSEOUT_PASS') && qualifications.failed_checks==0);
    assert(strcmp(qualifications.master_manifest_sha256,k1_pm_io('sha',fullfile(master,'manifest.json'))));
    assert(numel(qualifications.cases)==6);
end
paths={fullfile(master,'anchor')}; ids={'U6_N12_ANCHOR'}; winds=6; modes=12;
for j=1:6
    c=cfg.cases(j); root=fullfile(master,'cases',c.case_id); m=k1_pm_io('read',fullfile(root,'manifest.json')); fp=k1_pm_io('sha',fullfile(root,'manifest.json'));
    assert(strcmp(m.case_id,c.case_id) && m.seed==260001);
    k1_pm_io('verify_manifest',root); k1_pm_io('verify_seal',fullfile(root,'model_stage'),'committed.json',[fp ':MODELS:' fp]);
    pred=k1_pm_io('sha',fullfile(root,'model_stage','committed.json'));
    jobs={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
    for k=1:5
        dest=fullfile(root,'jobs',jobs{k}); k1_pm_io('verify_seal',dest,'committed.json',[pred ':C' num2str(k-1) ':' fp]);
        a=k1_pm_io('read',fullfile(dest,'audit.json')); g=k1_pm_io('read',fullfile(dest,'gates.json'));
        mark=k1_pm_io('read',fullfile(dest,'solve_invocation.json'));
        assert(a.failed_checks==0 && strcmp(a.status,'PASS') && strcmp(g.status,'PASS') && mark.invocations==1);
        pred=k1_pm_io('sha',fullfile(dest,'committed.json'));
    end
    if local_closeout
        qc=qualifications.cases(j); assert(strcmp(qc.case_id,c.case_id));
        pm_sensitivity_io('verify',qc.audit_root,qc.files);
        numeric=k1_pm_io('read',qc.numeric_receipt); graphics=k1_pm_io('read',qc.graphics_receipt);
        assert(numeric.failed_checks==0 && numeric.new_acoustic_solves==0 && ...
            contains(numeric.status,'PASS') && graphics.failed_checks==0 && graphics.checked_error_curves==12 && ...
            strcmp(graphics.status,'K1_SMALL_PANEL_FIG_AND_ZERO_PIXEL_PASS') && graphics.pixel_tolerance==0);
        original_work=fullfile(root,'finalization','attempts','attempt_0001','work');
        for field={'common.mat','comparison.mat'}
            assert(strcmp(k1_pm_io('sha',fullfile(qc.work,field{1})),k1_pm_io('sha',fullfile(original_work,field{1}))));
        end
        paths{end+1}=qc.work; %#ok<AGROW>
    else
        k1_pm_io('verify_seal',fullfile(root,'finalization'),'committed.json',[pred ':FINAL:' fp]);
        p=k1_pm_io('read',fullfile(root,'finalization','pointer.json')); cp=k1_pm_io('read',fullfile(root,'checkpoint.json'));
        assert(strcmp(cp.status,'K1_CASE_COMPLETE') && cp.failed_checks==0);
        paths{end+1}=k1_pm_io('path',fullfile(root,'finalization'),p.relative_work_dir); %#ok<AGROW>
    end
    ids{end+1}=c.case_id; winds(end+1)=c.U_PM_mps; modes(end+1)=c.modes; %#ok<AGROW>
end
metrics=table; points=table; surfaces=table; records=struct([]); inputs=struct('case_id',{},'comparison_sha256',{},'common_sha256',{});
for j=1:7
    q=load(fullfile(paths{j},'common.mat')); a=load(fullfile(paths{j},'comparison.mat'));
    audit=k1_pm_io('read',fullfile(paths{j},'audit','receipt.json'));
    final=k1_pm_io('read',fullfile(paths{j},'independent_final_metrics_audit.json'));
    assert(audit.failed_checks==0 && final.failed_checks==0 && a.reference_qualified && strcmp(a.reference_job,'bie_space'));
    assert(isequal(q.x,(-205:205)'*.1953125) && isequal(q.z,3*ones(411,1)) && q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000);
    assert(strcmp(q.cfg.phase_contract_version,'2.0.0') && all(a.bie_values<1e-4));
    b=a.U(:,2); scale=norm(b)/sqrt(411); mask=abs(b)>=.01*scale;
    assert(isequal(mask,a.mask) && abs(scale-a.scale)<1e-12);
    if j>1 && ~test_only
        expected_profile=readtable(fullfile(master,'cases',ids{j},'input','coefficients.csv'));
        assert(q.cfg.U_PM_mps==winds(j) && q.cfg.modes==modes(j) && strcmp(q.cfg.case_id,ids{j}) && isequal(q.profile,expected_profile));
    end
    item=struct('case_id',ids{j},'U_PM_mps',winds(j),'N',modes(j),'x',q.x,'U',a.U,'scale',scale,'mask',mask,'test_only',test_only);
    if isempty(records), records=item; else, records(j)=item; end
    M=readtable(fullfile(paths{j},'model_metrics.csv'));
    for model=[1 3]
        method="PS"; unc=a.punc; if model==3, method="FullFK"; unc=a.func; end
        u=a.U(:,model); residual=u-b; amplitude=abs(u)-abs(b); phase=angle(u.*conj(b)); shown=phase; shown(~mask)=NaN;
        E=norm(residual)/norm(b); threshold=5*norm(a.bunc+unc)/norm(b); ix=find(string(M.method)==method); assert(isscalar(ix));
        assert(abs(E-M.complex_relative_L2(ix))<1e-12 && abs(threshold-M.five_times_empirical_resolution(ix))<1e-12 && logical(M.resolved(ix))==(E>threshold));
        phase_energy=4*abs(u).*abs(b).*sin(phase/2).^2;
        assert(max(abs(abs(residual).^2-amplitude.^2-phase_energy))<1e-10);
        fraction=NaN; if sum(abs(residual).^2)>0, fraction=sum(phase_energy)/sum(abs(residual).^2); end
        [peak,index]=max(abs(residual));
        metrics=[metrics;table(string(ids{j}),260001,winds(j),modes(j),modes(j)*2*pi/160,method,E,threshold,E>threshold, ...
            norm(amplitude)/norm(b),fraction,nnz(mask),peak,q.x(index),norm(b-q.flat)/norm(q.flat), ...
            'VariableNames',{'case_id','seed','U_PM_mps','N','Kmax_rad_per_m','method','complex_relative_L2', ...
            'five_times_empirical_resolution','resolved','amplitude_relative_L2','phase_energy_fraction','valid_phase_points', ...
            'peak_complex_residual_Pa','peak_x_m','BIE_rough_flat_relative_L2'})]; %#ok<AGROW>
        points=[points;table(repmat(string(ids{j}),411,1),repmat(method,411,1),q.x,q.z,real(residual),imag(residual),abs(residual), ...
            abs(residual)/scale,amplitude,shown,mask,'VariableNames',{'case_id','method','x_m','z_m','residual_real_Pa','residual_imag_Pa', ...
            'residual_abs_Pa','residual_over_own_BIE_RMS','amplitude_difference_Pa','masked_wrapped_phase_rad','phase_valid'})]; %#ok<AGROW>
    end
    inputs(j)=struct('case_id',ids{j},'comparison_sha256',k1_pm_io('sha',fullfile(paths{j},'comparison.mat')), ...
        'common_sha256',k1_pm_io('sha',fullfile(paths{j},'common.mat')));
    if j==1, surface_file=fullfile(master,'anchor','K0_U6_N12_metrics.json');
    else, surface_file=fullfile(master,'cases',ids{j},'input','k0_metrics.json'); end
    s=k1_pm_io('read',surface_file);
    surface_row=table(string(ids{j}),winds(j),modes(j),s.Kmax_rad_per_m,s.dx_m,s.Knyquist_rad_per_m, ...
        s.raw_rms_height_m,s.embedded_rms_height_m,s.raw_rms_second_per_m,s.embedded_rms_second_per_m, ...
        s.variance_above_Kmax_m2,s.embedded_power_above_Kmax_m2,s.embedding_spectral_L1_power_change_m2, ...
        'VariableNames',{'case_id','U_PM_mps','N','physical_Kmax_rad_per_m','sampling_dx_m','sampling_Knyquist_rad_per_m', ...
        'raw_RMS_height_m','C2_RMS_height_m','raw_RMS_second_per_m','C2_RMS_second_per_m','PM_variance_above_Kmax_m2', ...
        'C2_power_above_raw_Kmax_m2','raw_C2_spectral_L1_power_difference_m2'});
    surface_row.continuous_height_variance_coverage=s.continuous_coverage_0_Kmax;
    surface_row.Kpeak_rad_per_m=s.Kpeak_rad_per_m;
    surface_row.raw_RMS_slope=s.raw_rms_slope;
    surface_row.C2_RMS_slope=s.embedded_rms_slope;
    surface_row.raw_max_abs_slope=s.raw_max_abs_slope;
    surface_row.C2_max_abs_slope=s.embedded_max_abs_slope;
    surface_row.statistics_samples=s.sample_count;
    surface_row.receiver_count=411;
    surface_row.receiver_dx_m=.1953125;
    surfaces=[surfaces;surface_row]; %#ok<AGROW>
end
metrics.BIE_origin=repmat("this_recovery_independent",height(metrics),1);
metrics.BIE_origin(metrics.case_id=="U6_N07")="imported_audited_N7";
metrics.BIE_origin(metrics.case_id=="U6_N12_ANCHOR")="qualified_R1_anchor";
writetable(metrics,fullfile(out,'seven_condition_metrics.csv')); writetable(points,fullfile(out,'point_errors.csv'));
writetable(surfaces,fullfile(out,'raw_C2_surface_spectral_differences.csv'));
contrasts=table;
for method=["PS" "FullFK"]
    t=metrics(metrics.method==method,:);
    low12=t.complex_relative_L2(t.U_PM_mps==5 & t.N==12); low40=t.complex_relative_L2(t.U_PM_mps==5 & t.N==40);
    high12=t.complex_relative_L2(t.U_PM_mps==8 & t.N==12); high40=t.complex_relative_L2(t.U_PM_mps==8 & t.N==40);
    interaction=(high40-high12)-(low40-low12);
    empirical_floor=sum(t.five_times_empirical_resolution(ismember(t.U_PM_mps,[5 8]) & ismember(t.N,[12 40])));
    contrasts=[contrasts;table(method,low40-low12,high40-high12,interaction,empirical_floor,abs(interaction)>empirical_floor, ...
        'VariableNames',{'method','fixed_U5_N40minusN12','fixed_U8_N40minusN12','interaction_U8minusU5_N40minusN12', ...
        'sum_four_five_times_empirical_floors','interaction_empirically_resolved'})]; %#ok<AGROW>
end
writetable(contrasts,fullfile(out,'paired_model_error_interactions.csv'));
paired=table;
pair_ids={'U6_N07','U6_N12_ANCHOR';'U6_N12_ANCHOR','U6_N40'; ...
    'U5_N12','U6_N12_ANCHOR';'U6_N12_ANCHOR','U8_N12';'U5_N12','U8_N12'; ...
    'U5_N40','U6_N40';'U6_N40','U8_N40';'U5_N40','U8_N40'; ...
    'U5_N12','U5_N40';'U8_N12','U8_N40'};
for method=["PS" "FullFK"]
    t=metrics(metrics.method==method,:);
    for pair=1:size(pair_ids,1)
        lo=t(t.case_id==string(pair_ids{pair,1}),:); hi=t(t.case_id==string(pair_ids{pair,2}),:);
        assert(height(lo)==1 && height(hi)==1);
        delta=hi.complex_relative_L2-lo.complex_relative_L2;
        floor=hi.five_times_empirical_resolution+lo.five_times_empirical_resolution;
        paired=[paired;table(method,lo.case_id,hi.case_id,delta,floor,abs(delta)>floor, ...
            'VariableNames',{'method','from_case','to_case','error_difference','sum_five_times_empirical_floors','difference_empirically_resolved'})]; %#ok<AGROW>
    end
end
writetable(paired,fullfile(out,'paired_model_error_contrasts.csv'));
field_changes=table;
for pair=1:size(pair_ids,1)
    lo=find(strcmp(ids,pair_ids{pair,1})); hi=find(strcmp(ids,pair_ids{pair,2}));
    for model=1:3
        names=["PS" "BIE" "FullFK"];
        change=norm(records(hi).U(:,model)-records(lo).U(:,model))/norm(records(lo).U(:,model));
        field_changes=[field_changes;table(string(pair_ids{pair,1}),string(pair_ids{pair,2}),names(model),change, ...
            'VariableNames',{'from_case','to_case','method','field_change_relative_to_from_field'})]; %#ok<AGROW>
    end
end
writetable(field_changes,fullfile(out,'paired_physical_field_changes.csv'));
k1_pm_io('write',fullfile(out,'input_versions.json'),inputs);
figures=render_k1_pm_summary_fields(records,fullfile(out,'figures')); assert(figures.failed_checks==0 && figures.test_only==test_only);
fid=fopen(fullfile(out,'report.md'),'w','n','UTF-8'); assert(fid>=0); guard=onCleanup(@()fclose(fid));
if test_only
    fprintf(fid,'# TEST_ONLY 汇总软件夹具\n\n仅重放历史R1数组验证汇总程序，没有六个新K1声学结果。所有模拟检查点只属于测试夹具，不是数值资格。\n\n');
else
    fprintf(fid,'# K1 风速—物理谱截止联合声场敏感性\n\nN7从已审查旧K1任务导入；其余五面在本恢复任务独立完成原PS/BIE/FullFK层级。U6/N12为严格匹配的历史R1锚点。\n\n');
end
fprintf(fid,'固定0.10m二维条带、4kHz、seed260001、原C2、phase contract2.0.0、411真实接收点。每个海面使用自己的BIE范数与0.01全线RMS相位阈值。N40不表示完整PM谱已收敛。\n\n');
fprintf(fid,'| U_PM (m/s) | N | 模型 | 复场L2 (%%) | 5倍经验阈值 (%%) | 可分辨 |\n|---:|---:|---|---:|---:|---|\n');
for j=1:height(metrics), fprintf(fid,'|%g|%d|%s|%.8g|%.8g|%d|\n',metrics.U_PM_mps(j),metrics.N(j),metrics.method(j),100*metrics.complex_relative_L2(j),100*metrics.five_times_empirical_resolution(j),metrics.resolved(j)); end
fprintf(fid,'\n可分辨性是经验数值分辨，不是严格真解或工程合格线。未分辨比较如实保留，不自动调网格/门槛。七个条件为同一随机块，411点不当独立重复，不外推总体显著性、宽带/CIR或通信。原始/C2几何与谱差异保存在每案例input/k0_*文件，不把它们当相同带限面。\n');
for model={'PS','FullFK'}
    for kind={'absolute','normalized','amplitude','phase'}
        fprintf(fid,'\n![%s %s](figures/K1_%s_%s.png)\n',model{1},kind{1},model{1},kind{1});
    end
end
if local_closeout
    fprintf(fid,'\n资格由原MODELS/C0–C4封存链、本地独立保存数值审计及新固定画布图件零像素审计组成。旧大画布图件FAIL及原远端退出回执保留，未伪造旧FINAL封存；见local-closeout资格清单。\n');
end
clear guard
receipt=struct('status','K1_SUMMARY_AUDITED_PASS','failed_checks',0,'new_cases',6,'conditions_including_anchor',7, ...
    'model_rows',height(metrics),'point_rows',height(points),'unresolved_comparisons',nnz(~metrics.resolved), ...
    'new_BIE_invocations',25,'total_BIE_invocations_including_import',30,'imported_BIE_invocations',5,'figure_curves',figures.curves_checked,'elapsed_s',toc(timer), ...
    'scope','single seed, single frequency, same 411 actual receivers; no wideband or communications');
receipt.test_only=test_only;
receipt.qualification_mode=mode;
if test_only, receipt.status='K1_SUMMARY_TEST_ONLY_PASS'; receipt.new_BIE_invocations=0; receipt.total_BIE_invocations_including_import=0; receipt.imported_BIE_invocations=0; receipt.new_physical_cases=0; end
k1_pm_io('write',fullfile(out,'receipt.json'),receipt); disp(jsonencode(receipt));
end
