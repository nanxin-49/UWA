function receipt=summarize_r2_pm_fields(work_dirs,output_dir)
%SUMMARIZE_R2_PM_FIELDS Saved-data analysis for one seed or ordered R1+R2.
% No solver, no statistics treating receiver points as independent samples.
arguments
    work_dirs (1,:) cell
    output_dir (1,:) char
end
assert(ismember(numel(work_dirs),[1 3]) && ~isfolder(output_dir));
summary=table; points=table; regions=table; norms=table; inputs=struct('path',{},'sha256',{});
anchor=[]; seeds=zeros(1,numel(work_dirs));
for j=1:numel(work_dirs)
    root=work_dirs{j}; q=load(fullfile(root,'common.mat')); a=load(fullfile(root,'comparison.mat'));
    audit=r2_pm_io('read',fullfile(root,'independent_final_metrics_audit.json'));
    kernel=r2_pm_io('read',fullfile(root,'audit','receipt.json'));
    assert(strcmp(audit.status,'PASS') && audit.failed_checks==0 && kernel.failed_checks==0 && a.reference_qualified);
    assert(strcmp(a.reference_job,'bie_space') && all(a.bie_values<1e-4));
    seeds(j)=q.cfg.seed; assert(ismember(seeds(j),[260001 260002 260003]));
    if isempty(anchor), anchor=q;
    else
        assert(isequal(q.x,anchor.x) && isequal(q.z,anchor.z) && isequal(q.cfg.source,anchor.cfg.source));
        assert(isequal(q.profile.k_rad_per_m,anchor.profile.k_rad_per_m) && isequal(q.profile.Sk_m3,anchor.profile.Sk_m3));
        assert(isequal(q.cfg.ps_recipe,anchor.cfg.ps_recipe) && isequal(q.cfg.fk_recipe,anchor.cfg.fk_recipe));
    end
    assert(numel(q.x)==411 && all(q.z==3) && q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000);
    assert(strcmp(q.cfg.phase_contract_version,'2.0.0') && q.cfg.phase_mask_fraction_of_BIE_RMS==.01);
    T=readtable(fullfile(root,'receiver_fields_and_errors.csv')); M=readtable(fullfile(root,'model_metrics.csv'));
    assert(max(abs(reshape(T.real_Pa+1i*T.imag_Pa,411,3)-a.U),[],'all')<1e-11);
    for f={'common.mat','comparison.mat','receiver_fields_and_errors.csv','model_metrics.csv', ...
            'independent_final_metrics_audit.json','audit/receipt.json'}
        path=fullfile(root,f{1}); inputs(end+1)=struct('path',path,'sha256',r2_pm_io('sha',path)); %#ok<AGROW>
    end
    b=a.U(:,2); scale=norm(b)/sqrt(411); mask=abs(b)>=.01*scale;
    assert(isequal(mask,a.mask) && abs(scale-a.scale)<1e-12);
    for m=[1 3]
        method="PS"; unc=a.punc; if m==3, method="FullFK"; unc=a.func; end
        u=a.U(:,m); e=u-b; amp=abs(u)-abs(b); phase=angle(u.*conj(b)); sq=abs(e).^2;
        phase_energy=4*abs(u).*abs(b).*sin(phase/2).^2;
        assert(max(abs(sq-amp.^2-phase_energy))<1e-10);
        E=norm(e)/norm(b); limit=5*norm(a.bunc+unc)/norm(b);
        row=find(string(M.method)==method); assert(isscalar(row));
        assert(abs(E-M.complex_relative_L2(row))<1e-12 && abs(limit-M.five_times_empirical_resolution(row))<1e-12);
        [peak,ix]=max(abs(e)); [amp_peak,ia]=max(abs(amp));
        valid=find(mask); maxphase=NaN; ip=NaN; phase_rms=NaN;
        if ~isempty(valid)
            [maxphase,jp]=max(abs(phase(valid))); ip=q.x(valid(jp));
            phase_rms=sqrt(sum(abs(b(mask)).^2.*phase(mask).^2)/sum(abs(b(mask)).^2));
        end
        energy=sum(sq); fraction=NaN; topfraction=NaN;
        [~,order]=sort(sq,'descend');
        if energy>0, fraction=sum(phase_energy)/energy; topfraction=sum(sq(order(1:41)))/energy; end
        summary=[summary;table(seeds(j),method,E,limit,E>limit,norm(amp)/norm(b), ...
            fraction,phase_rms,maxphase,ip,nnz(mask),peak,q.x(ix),amp_peak,q.x(ia),topfraction, ...
            norm(b-q.flat)/norm(q.flat),'VariableNames',{'seed','method','complex_relative_L2', ...
            'five_times_empirical_resolution','resolved','amplitude_relative_L2','phase_energy_fraction', ...
            'BIE_power_weighted_phase_RMS_rad','max_abs_phase_rad','max_phase_x_m','valid_phase_points', ...
            'peak_complex_residual_Pa','peak_complex_x_m','max_abs_amplitude_difference_Pa','max_amplitude_x_m', ...
            'top_41_energy_fraction','BIE_rough_flat_relative_L2'})]; %#ok<AGROW>
        shown_phase=phase; shown_phase(~mask)=NaN;
        points=[points;table(repmat(seeds(j),411,1),repmat(method,411,1),q.x,q.z, ...
            real(e),imag(e),abs(e),abs(e)/scale,amp,shown_phase,mask,amp.^2,phase_energy, ...
            'VariableNames',{'seed','method','x_m','z_m','residual_real_Pa','residual_imag_Pa','residual_abs_Pa', ...
            'residual_over_BIE_RMS','amplitude_difference_Pa','phase_difference_rad','phase_valid', ...
            'amplitude_error_energy_Pa2','phase_error_energy_Pa2'})]; %#ok<AGROW>
        bins=[-Inf -20 0 20 Inf];
        for k=1:4
            ii=q.x>=bins(k) & q.x<bins(k+1); share=NaN;
            if energy>0, share=sum(sq(ii))/energy; end
            regions=[regions;table(seeds(j),method,k,min(q.x(ii)),max(q.x(ii)),nnz(ii), ...
                sum(sq(ii)),share,sum(amp(ii).^2),sum(phase_energy(ii)), ...
                'VariableNames',{'seed','method','region','left_x_m','right_x_m','points', ...
                'residual_energy_Pa2','residual_energy_fraction','amplitude_energy_Pa2','phase_energy_Pa2'})]; %#ok<AGROW>
        end
        edges=round(linspace(1,411,5));
        for k=1:4
            ii=edges(k):edges(k+1);
            norms=[norms;table(seeds(j),method,k,q.x(ii(1)),q.x(ii(end)),norm(e(ii))/norm(b(ii)), ...
                'VariableNames',{'seed','method','region','left_x_m','right_x_m','complex_relative_L2'})]; %#ok<AGROW>
        end
    end
end
if numel(work_dirs)==3, assert(isequal(seeds,[260001 260002 260003])); end
mkdir(output_dir);
writetable(summary,fullfile(output_dir,'seed_metrics.csv')); writetable(points,fullfile(output_dir,'point_errors.csv'));
writetable(regions,fullfile(output_dir,'disjoint_energy_regions.csv')); writetable(norms,fullfile(output_dir,'shared_endpoint_region_norms.csv'));
ordering=table;
for seed=seeds
    ps=summary(summary.seed==seed & summary.method=="PS",:); fk=summary(summary.seed==seed & summary.method=="FullFK",:);
    ordering=[ordering;table(seed,ps.complex_relative_L2/fk.complex_relative_L2, ...
        ps.complex_relative_L2>fk.complex_relative_L2,ps.resolved,fk.resolved, ...
        'VariableNames',{'seed','PS_over_FullFK_L2','FullFK_smaller','PS_resolved','FullFK_resolved'})]; %#ok<AGROW>
end
writetable(ordering,fullfile(output_dir,'model_ordering.csv'));
receipt=struct('status','SAVED_DATA_SUMMARY_PASS','failed_checks',0,'seeds',seeds,'inputs',inputs, ...
    'new_acoustic_solves',0,'phase_mask_fraction_of_seed_BIE_RMS',.01, ...
    'interpretation','Three predetermined surfaces only; no confidence interval, no engineering error threshold, no iid receiver inference');
r2_pm_io('write',fullfile(output_dir,'receipt.json'),receipt);
fid=fopen(fullfile(output_dir,'report.md'),'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'# R2 saved-field comparison\n\nFixed 0.10 m 2-D strip, 4 kHz, 411 physical receivers.\n\n');
fprintf(fid,'| Seed | Model | Complex L2 | Five-times resolution | Phase energy fraction | Peak x (m) |\n|---|---|---:|---:|---:|---:|\n');
for k=1:height(summary)
    fprintf(fid,'| %d | %s | %.12g | %.12g | %.9g | %.9g |\n',summary.seed(k),summary.method(k), ...
        summary.complex_relative_L2(k),summary.five_times_empirical_resolution(k),summary.phase_energy_fraction(k),summary.peak_complex_x_m(k));
end
fprintf(fid,'\nEnergy regions are disjoint; convergence/norm regions retain R1 shared endpoints. Phase decomposition uses all points; displayed phase and power-weighted RMS use the fixed mask.\n');
fprintf(fid,'\nRead model_ordering.csv per seed; an unresolved model difference does not establish equivalence. These three fixed surfaces do not establish a population law.\n');
end
