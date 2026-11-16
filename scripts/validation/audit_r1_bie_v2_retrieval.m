function receipt=audit_r1_bie_v2_retrieval(source,out)
%AUDIT_R1_BIE_V2_RETRIEVAL Saved-data review; no density/PS/FullFK solves.
arguments
    source (1,:) char
    out (1,:) char
end
source=char(java.io.File(source).getCanonicalPath()); out=char(java.io.File(out).getCanonicalPath());
assert(~startsWith(out,[source filesep]) && ~strcmp(out,source));
assert(isfile(fullfile(out,'retrieval_integrity.json')) && ~isfolder(fullfile(out,'reaudit')));
previous_path=path; path_guard=onCleanup(@()path(previous_path)); addpath(fullfile(source,'code'),'-begin');
threads=maxNumCompThreads(1); thread_guard=onCleanup(@()maxNumCompThreads(threads));
timer=tic; diary(fullfile(out,'diary.log')); diary_guard=onCleanup(@()diary('off'));
try
    integrity=r1_bie_continuation_io('read',fullfile(out,'retrieval_integrity.json'));
    assert(strcmp(integrity.status,'PASS') && integrity.snapshots_task_files_byte_identical);
    chain=r1_bie_continuation_io('verify_migration',source);
    c4=fullfile(source,'jobs','bie_window');
    sig=[chain.predecessor ':C4:' chain.manifest_sha256 ':' chain.policy_sha256];
    assert(strcmp(r1_bie_continuation_io('c4_state',c4,sig),'SKIP_AUDITED'));
    final=fullfile(source,'finalization'); pointer=r1_bie_continuation_io('read',fullfile(final,'pointer.json'));
    work=r1_bie_continuation_io('path',final,pointer.relative_work_dir);
    fsig=[r1_bie_continuation_io('sha',fullfile(c4,'committed.json')) ':FINAL:' chain.manifest_sha256 ':' chain.policy_sha256];
    seal=r1_bie_continuation_io('verify_seal',final,'committed.json',fsig,{'pointer.json'});
    assert(strcmp(seal.status,'AUDITED_PASS') && strcmp(pointer.status,'PASS'));
    q=load(fullfile(work,'common.mat')); a=load(fullfile(work,'comparison.mat'));
    assert(strcmp(r1_bie_continuation_io('sha',fullfile(work,'common.mat')),r1_bie_continuation_io('sha',fullfile(source,'common.mat'))));
    assert(q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000 && ...
        strcmp(q.cfg.phase_contract_version,'2.0.0') && numel(q.x)==411 && all(q.z==3));
    assert(a.reference_qualified && strcmp(a.reference_job,'bie_space') && isequal(size(a.U),[411 3]));
    ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'}; B=cell(5,1);
    for j=1:5
        path1=fullfile(source,'jobs',ids{j},'data.mat'); path2=fullfile(work,'jobs',ids{j},'data.mat');
        assert(strcmp(r1_bie_continuation_io('sha',path1),r1_bie_continuation_io('sha',path2)));
        B{j}=load(path1);
    end
    assert(isequal(a.U(:,2),B{4}.field));
    % Every solve entrypoint is shadowed while saved-density/kernel audits run.
    guards=fullfile(out,'solve_sentinels'); mkdir(guards);
    for name={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
            'run_pe_1d_surface_reflection_validation','closeout_s3a_pm_fullfk'}
        fid=fopen(fullfile(guards,[name{1} '.m']),'w'); assert(fid>=0);
        fprintf(fid,'function varargout=%s(varargin)\nerror(''r1review:ForbiddenSolve'',''Only saved-data audits permitted'');\nend\n',name{1}); fclose(fid);
    end
    addpath(guards,'-begin'); rehash;
    assert(strcmp(which('solve_helmholtz_bie_halfplane_vertical'),fullfile(guards,'solve_helmholtz_bie_halfplane_vertical.m')));
    review=fullfile(out,'reaudit'); mkdir(review);
    for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv'}
        copyfile(fullfile(work,name{1}),fullfile(review,name{1}));
    end
    copyfile(fullfile(work,'input'),fullfile(review,'input'));
    for j=1:5
        dest=fullfile(review,'jobs',ids{j}); mkdir(dest);
        copyfile(fullfile(work,'jobs',ids{j},'data.mat'),fullfile(dest,'data.mat'));
    end
    c4_review=fullfile(review,'jobs','bie_window'); copyfile(fullfile(c4,'points.csv'),fullfile(c4_review,'points.csv'));
    c4_audit=audit_r1_bie_manual_stage(c4_review,fullfile(review,'common.mat'),[70 12 96]);
    all_audit=audit_piston_pm_r1(review); diary(fullfile(out,'diary.log'));
    final_audit=audit_r1_bie_continuation_final(review,source);
    assert(c4_audit.failed_checks==0 && all_audit.failed_checks==0 && final_audit.failed_checks==0);
    % Fresh FIG data checks; compare the two already-exported PNGs per pair,
    % avoiding cross-platform font/render changes or rewriting original figures.
    figures=fullfile(work,'figures_final'); checked=0;
    for name={'F02_pointwise_errors','F04_FullFK_pointwise_errors'}
        fig=openfig(fullfile(figures,[name{1} '.fig']),'new','invisible'); fig_guard=onCleanup(@()close(fig));
        lines=findall(fig,'Type','line');
        for j=1:numel(lines)
            info=lines(j).UserData;
            if ~isstruct(info) || ~isfield(info,'panel'), continue; end
            assert(ismember(info.method,{'PS','FullFK'}));
            m=3; if strcmp(info.method,'PS'), m=1; end
            e=a.U(:,m)-a.U(:,2);
            switch info.panel
                case 1, expected=abs(e);
                case 2, expected=abs(e)/a.scale;
                case 3, expected=abs(a.U(:,m))-abs(a.U(:,2));
                case 4, expected=angle(a.U(:,m).*conj(a.U(:,2))); expected(~a.mask)=NaN;
                otherwise, error('Unexpected saved figure panel');
            end
            actual=lines(j).YData(:);
            assert(isequal(lines(j).XData(:),q.x) && isequal(isnan(actual),isnan(expected)));
            ii=~isnan(expected); assert(max(abs(actual(ii)-expected(ii)))<1e-11); checked=checked+1;
        end
        clear fig_guard;
    end
    assert(checked==12);
    im1=imread(fullfile(figures,'F02_pointwise_errors.png')); im2=imread(fullfile(figures,'F02_audit_reexport.png'));
    assert(isequal(im1,im2));
    im1=imread(fullfile(figures,'F04_FullFK_pointwise_errors.png')); im2=imread(fullfile(figures,'F04_audit_reexport.png'));
    assert(isequal(im1,im2));
    % Independent along-line metrics and a stable amplitude/phase energy identity.
    T=readtable(fullfile(work,'model_metrics.csv')); metrics=table; regions=table; errors=table;
    edges=round(linspace(1,411,5));
    for m=[1 3]
        method="PS"; model_unc=a.punc; if m==3, method="FullFK"; model_unc=a.func; end
        ref=a.U(:,2); u=a.U(:,m); e=u-ref; amplitude=abs(u)-abs(ref); phase=angle(u.*conj(ref));
        E=norm(e)/norm(ref); resolution=5*norm(a.bunc+model_unc)/norm(ref);
        row=find(string(T.method)==method); assert(abs(E-T.complex_relative_L2(row))<1e-12);
        amp_energy=sum(amplitude.^2); phase_energy=sum(4*abs(u).*abs(ref).*sin(phase/2).^2);
        energy=sum(abs(e).^2); assert(abs(amp_energy+phase_energy-energy)/energy<1e-10);
        [peak,ix]=max(abs(e)); [peak_amp,ia]=max(abs(amplitude)); [peak_phase,ip]=max(abs(phase));
        phase_rms=sqrt(sum(abs(ref).^2.*phase.^2)/sum(abs(ref).^2));
        sorted=sort(abs(e).^2,'descend');
        metrics=[metrics;table(method,E,resolution,E/resolution,peak,q.x(ix),norm(amplitude)/norm(ref), ...
            peak_amp,q.x(ia),phase_rms,peak_phase,q.x(ip),phase_energy/energy, ...
            sum(abs(e(q.x>=0)).^2)/energy,sum(sorted(1:41))/energy, ...
            'VariableNames',{'method','complex_relative_L2','five_times_resolution','error_to_resolution_ratio', ...
            'peak_residual_Pa','peak_x_m','amplitude_relative_L2','peak_amplitude_difference_Pa','peak_amplitude_x_m', ...
            'reference_weighted_phase_RMS_rad','peak_phase_abs_rad','peak_phase_x_m','phase_energy_fraction', ...
            'nonnegative_x_energy_fraction','largest_41_points_energy_fraction'})]; %#ok<AGROW>
        for j=1:4
            ii=edges(j):edges(j+1); Epart=norm(e(ii))/norm(ref(ii));
            unique_ii=edges(j):edges(j+1)-1; if j==4, unique_ii=edges(j):edges(j+1); end
            regions=[regions;table(method,j,q.x(ii(1)),q.x(ii(end)),Epart,sum(abs(e(unique_ii)).^2)/energy, ...
                'VariableNames',{'method','quarter','left_x_m','right_x_m','relative_L2_shared_endpoints','energy_fraction_disjoint_bins'})]; %#ok<AGROW>
        end
        errors=[errors;table(repmat(method,411,1),q.x,q.z,real(e),imag(e),abs(e),amplitude,phase, ...
            5*(a.bunc+model_unc),'VariableNames',{'method','x_m','z_m','residual_real_Pa','residual_imag_Pa', ...
            'residual_abs_Pa','amplitude_difference_Pa','phase_difference_rad','five_times_point_empirical_budget_Pa'})]; %#ok<AGROW>
    end
    % Verify historical pairing with native saved fields and actual coordinates.
    gaussian_ids={'main96','self64','space14','window60','window70','pe','fullfk','pe_grid'}; G=cell(8,1);
    for j=1:8
        g=load(fullfile(source,'input','gaussian',[gaussian_ids{j} '.mat']),'payload'); G{j}=g.payload;
        assert(isequal(G{j}.x_m(1:411),q.x) && isequal(G{j}.z_m(1:411),q.z));
    end
    profile=readtable(fullfile(source,'input','seed_260001_coefficients.csv')); assert(isequal(profile,q.profile));
    gu=[G{6}.u_native(1:411) G{1}.u_native(1:411) G{7}.u_native(1:411)]; assert(isequal(gu,a.gaussian.fields));
    gaussian_E=[norm(gu(:,1)-gu(:,2));norm(gu(:,3)-gu(:,2))]/norm(gu(:,2));
    paired=table(metrics.method,metrics.complex_relative_L2,gaussian_E,metrics.complex_relative_L2./gaussian_E, ...
        'VariableNames',{'method','piston_own_BIE_relative_L2','Gaussian_own_BIE_relative_L2','conditional_error_ratio'});
    % Check why the exceptionally small window difference is a genuine L70 result.
    bn=B{5}.boundary; bn1=B{2}.boundary;
    window=struct('C1_L_m',bn1.half_width_m,'C4_L_m',bn.half_width_m,'C1_nodes',bn1.node_count,'C4_nodes',bn.node_count, ...
        'C1_panel_width_m',2*bn1.half_width_m/bn1.panel_count,'C4_panel_width_m',2*bn.half_width_m/bn.panel_count, ...
        'C4_x_min_m',min(bn.node_x_m),'C4_x_max_m',max(bn.node_x_m), ...
        'C4_density_outside_50m_fraction',norm(bn.density(abs(bn.node_x_m)>=50))/norm(bn.density), ...
        'C4_density_outside_60m_fraction',norm(bn.density(abs(bn.node_x_m)>=60))/norm(bn.density), ...
        'data_files_distinct',~strcmp(r1_bie_continuation_io('sha',fullfile(source,'jobs','bie_main','data.mat')), ...
        r1_bie_continuation_io('sha',fullfile(c4,'data.mat'))));
    assert(window.data_files_distinct && bn.half_width_m==70 && bn.node_count==4480 && bn1.node_count==3840);
    R=readtable(fullfile(source,'sessions','20261116T091043_470983242_3891066','resources.csv'), ...
        'ReadVariableNames',false,'Delimiter',',');
    available=R.Var2; rss=R.Var3;
    runtime=r1_bie_continuation_io('read',fullfile(source,'sessions','20261116T091043_470983242_3891066','runtime.json'));
    resources=struct('minimum_available_bytes',min(available),'baseline_bytes',runtime.RAM_baseline_bytes, ...
        'minimum_baseline_fraction',min(available)/runtime.RAM_baseline_bytes,'peak_sampled_group_RSS_KiB',max(rss));
    assert(resources.minimum_baseline_fraction>=.3);
    writetable(metrics,fullfile(out,'review_model_metrics.csv')); writetable(regions,fullfile(out,'review_regions.csv'));
    writetable(errors,fullfile(out,'review_point_errors.csv')); writetable(paired,fullfile(out,'conditional_Gaussian_pair.csv'));
    save(fullfile(out,'review_analysis.mat'),'metrics','regions','errors','paired','window','resources','-v7');
    receipt=struct('status','SAVED_RESULT_REVIEW_PASS','failed_checks',0,'elapsed_s',toc(timer), ...
        'source',source,'final_work',work,'reference_stage',3,'reference_scope','fixed 10cm 4kHz seed260001 empirical hierarchy', ...
        'c4_audit',c4_audit,'comprehensive_audit',all_audit,'reference_resolution_audit',final_audit, ...
        'checked_FIG_error_curves',checked,'saved_PNG_reexport_pixel_differences',[0 0], ...
        'metrics',table2struct(metrics),'conditional_pair',table2struct(paired),'window_diagnostics',window, ...
        'resources',resources,'max_threads',maxNumCompThreads,'new_acoustic_solves',0,'SSH_connections',0, ...
        'legacy_nonmonotonic_FAIL_preserved',true,'R2_authorized',false);
    r1_bie_continuation_io('write',fullfile(out,'review_receipt.json'),receipt);
    diary('off'); clear diary_guard;
    disp(jsonencode(struct('status',receipt.status,'failed_checks',0,'elapsed_s',receipt.elapsed_s, ...
        'metrics',receipt.metrics,'window',window,'new_acoustic_solves',0)));
catch ME
    r1_bie_continuation_io('write',fullfile(out,'failure.json'),struct('status','FAIL','identifier',ME.identifier,'error',ME.message));
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
