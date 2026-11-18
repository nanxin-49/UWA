function receipt = review_pm_wind_cutoff_saved(source_dir, output_dir)
%REVIEW_PM_WIND_CUTOFF_SAVED Independent local audit of retrieved K0, no solve.
% Copies only four per-case numerical files; never writes to the retrieval.
arguments
    source_dir (1,:) char
    output_dir (1,:) char
end
pm_sensitivity_io('absolute',source_dir); pm_sensitivity_io('absolute',output_dir);
assert(~isfolder(output_dir)); mkdir(output_dir);
logfile=fullfile(output_dir,'review.log'); timer=tic; cpu0=cputime;
threads=maxNumCompThreads(1); tg=onCleanup(@()maxNumCompThreads(threads));
[~,system_info]=memory; baseline=system_info.PhysicalMemory.Available;
transfer=jsondecode(fileread(fullfile(source_dir,'retrieval','manifest.json')));
pm_sensitivity_io('verify',source_dir,transfer);
binding=pm_sensitivity_io('sha',fullfile(source_dir,'manifest.json'));
assert(strcmp(binding,'cf48069437260281d34a560ba6c3bcee1047d61acdc18154e8c1d68aeb6fd6c6'));
% Use exactly the frozen independent auditor and IO helper.
assert(strcmp(pm_sensitivity_io('sha',which('audit_pm_wind_cutoff_k0_case')), ...
    pm_sensitivity_io('sha',fullfile(source_dir,'code','audit_pm_wind_cutoff_k0_case.m'))));
assert(strcmp(pm_sensitivity_io('sha',which('pm_sensitivity_io')), ...
    pm_sensitivity_io('sha',fullfile(source_dir,'code','pm_sensitivity_io.m'))));
cfg=jsondecode(fileread(fullfile(source_dir,'config.json')));
stats=readtable(fullfile(source_dir,'output','summary_attempt_001','statistics.csv'));
assert(height(stats)==45 && strcmp(cfg.kind,'K0_SURFACE_ONLY'));
audits=table; bounds=table; minimum=baseline; max_metric_difference=0;
try
    for j=1:numel(cfg.cases)
        [~,mem]=memory; minimum=min(minimum,mem.PhysicalMemory.Available);
        assert(minimum>=.3*baseline,'Local RAM reserve stop');
        assert(toc(timer)<180,'Local saved-review time stop');
        c=cfg.cases(j); id=sprintf('seed%d_U%d_N%02d',c.seed,c.U_PM_mps,c.N);
        src=fullfile(source_dir,'output',id); dest=fullfile(output_dir,id); mkdir(dest);
        commit=jsondecode(fileread(fullfile(src,'checkpoint.json')));
        assert(strcmp(commit.status,'AUDITED_PASS') && strcmp(commit.manifest_sha256,binding));
        pm_sensitivity_io('verify',src,commit.artifacts_sha256);
        for file={'surface.mat','surface.csv','coefficients.csv','spectra.csv'}
            copyfile(fullfile(src,file{1}),fullfile(dest,file{1}));
            assert(strcmp(pm_sensitivity_io('sha',fullfile(src,file{1})),pm_sensitivity_io('sha',fullfile(dest,file{1}))));
        end
        t=tic; a=audit_pm_wind_cutoff_k0_case(dest); elapsed=toc(t);
        q=load(fullfile(dest,'surface.mat'),'data'); d=q.data; s=d.samples; coef=d.coefficients;
        % Audit metric fields omitted by the frozen core auditor.
        actual=[mean(s.raw_elevation_m),min(s.raw_elevation_m),max(s.raw_elevation_m), ...
            mean(s.embedded_elevation_m),min(s.embedded_elevation_m),max(s.embedded_elevation_m), ...
            max(abs(s.raw_second_per_m)),max(abs(s.embedded_second_per_m)), ...
            sqrt(mean((s.raw_second_per_m./(1+s.raw_slope.^2).^1.5).^2)), ...
            sqrt(mean((s.embedded_second_per_m./(1+s.embedded_slope.^2).^1.5).^2))];
        names={'raw_mean_height_m','raw_min_elevation_m','raw_max_elevation_m', ...
            'embedded_mean_height_m','embedded_min_elevation_m','embedded_max_elevation_m', ...
            'raw_max_abs_second_per_m','embedded_max_abs_second_per_m', ...
            'raw_rms_geometric_curvature_per_m','embedded_rms_geometric_curvature_per_m'};
        for k=1:numel(names), assert(abs(actual(k)-d.metrics.(names{k}))<1e-12); end
        ix=stats.seed==c.seed & stats.U_PM_mps==c.U_PM_mps & stats.N==c.N; assert(nnz(ix)==1);
        fields=fieldnames(d.metrics);
        for k=1:numel(fields)
            delta=abs(stats.(fields{k})(ix)-d.metrics.(fields{k}));
            assert(delta<1e-11); max_metric_difference=max(max_metric_difference,delta);
        end
        % A certified continuous bound, not a sampled-only clearance claim.
        % |eta_raw| <= B0; |eta_raw'| <= B1; max|chi'|=1.875/8.
        amplitudes=hypot(coef.cos_depth_coeff_m,coef.sin_depth_coeff_m);
        B0=sum(amplitudes); B1=sum(coef.k_rad_per_m.*amplitudes);
        slope_bound=B1+(1.875/8)*B0;
        margin=(160/8192)/2*slope_bound;
        low=min(s.embedded_elevation_m)-margin; high=max(s.embedded_elevation_m)+margin;
        bounds=[bounds;table(c.seed,c.U_PM_mps,c.N,B0,slope_bound,margin,low,high, ...
            2-high,3+low,100+low,'VariableNames',{'seed','U_PM_mps','N','Fourier_amplitude_bound_m', ...
            'C2_slope_absolute_bound','grid_to_continuum_margin_m','C2_elevation_lower_bound_m', ...
            'C2_elevation_upper_bound_m','certified_auxiliary_plane_gap_m', ...
            'certified_receiver_clearance_m','certified_source_depth_clearance_m'})]; %#ok<AGROW>
        audits=[audits;table(string(id),a.reconstruction_max_abs,max(a.grid_rms_relative), ...
            a.grid_max_slope_absolute,elapsed,'VariableNames',{'case_id','reconstruction_max_abs', ...
            'grid_RMS_relative_max','grid_max_slope_absolute','elapsed_s'})]; %#ok<AGROW>
        pm_sensitivity_io('log',logfile,sprintf('%02d/45 %s local_saved_audit=PASS elapsed_s=%.4f',j,id,elapsed));
    end
    writetable(audits,fullfile(output_dir,'local_audit_summary.csv'));
    writetable(bounds,fullfile(output_dir,'continuous_geometry_bounds.csv'));
    figure_receipt=local_figures(source_dir,stats,output_dir);
    pm_sensitivity_io('verify',source_dir,transfer);
    receipt=struct('status','K0_LOCAL_SAVED_REVIEW_PASS','case_count',height(audits), ...
        'failed_checks',0,'matlab_release',version('-release'),'execution_mode','existing_local_MATLAB_MCP', ...
        'elapsed_s',toc(timer),'cpu_s',cputime-cpu0,'manifest_sha256',binding, ...
        'source_files_SHA_before_after',numel(transfer),'local_reconstruction_max_abs',max(audits.reconstruction_max_abs), ...
        'local_grid_RMS_relative_max',max(audits.grid_RMS_relative_max), ...
        'local_grid_max_slope_absolute',max(audits.grid_max_slope_absolute), ...
        'summary_metric_max_abs_difference',max_metric_difference, ...
        'minimum_available_RAM_bytes',minimum,'starting_available_RAM_bytes',baseline, ...
        'original_FIG_curve_count',figure_receipt.curves,'FIG_audit_path','figure_audit.json', ...
        'minimum_certified_auxiliary_plane_gap_m',min(bounds.certified_auxiliary_plane_gap_m), ...
        'minimum_certified_receiver_clearance_m',min(bounds.certified_receiver_clearance_m), ...
        'new_surface_generation_calls',0,'new_acoustic_solves',0,'SSH_connections',0, ...
        'next_allowed_step','K0 conclusions and K1 design; no K1 execution');
    pm_sensitivity_io('write',fullfile(output_dir,'receipt.json'),receipt);
    disp(jsonencode(receipt));
catch err
    pm_sensitivity_io('log',logfile,getReport(err,'extended'));
    pm_sensitivity_io('write',fullfile(output_dir,'failure.json'),struct('status','FAILED_REVIEW','error',err.message,'elapsed_s',toc(timer)));
    rethrow(err)
end
end

function receipt=local_figures(source_dir,stats,output_dir)
summary=fullfile(source_dir,'output','summary_attempt_001');
fields={'raw_rms_height_m','raw_rms_slope','raw_rms_second_per_m','embedded_rms_height_m', ...
    'embedded_rms_slope','embedded_rms_second_per_m','continuous_coverage_0_Kmax', ...
    'variance_above_Kmax_m2','embedded_power_above_Kmax_m2'};
expected=struct('x',{},'y',{}); n=0;
for k=1:numel(fields)
    for u=[5 6 8]
        for seed=260001:260003
            rows=stats(stats.U_PM_mps==u & stats.seed==seed,:); rows=sortrows(rows,'N');
            n=n+1; expected(n)=struct('x',rows.Kmax_rad_per_m(:),'y',rows.(fields{k})(:));
        end
    end
end
n1=local_match_fig(fullfile(summary,'statistics.fig'),expected);
expected=struct('x',{},'y',{}); n=0;
for u=[5 6 8]
    for N=[7 12 40]
        folder=fullfile(source_dir,'output',sprintf('seed260001_U%d_N%02d',u,N));
        s=readtable(fullfile(folder,'surface.csv')); p=readtable(fullfile(folder,'spectra.csv'));
        for e=0:1
            if e==0, y=s.raw_elevation_m; power=p.raw_bin_power_m2;
            else, y=s.embedded_elevation_m; power=p.embedded_bin_power_m2; end
            n=n+1; expected(n)=struct('x',s.x_m(:),'y',y(:));
            keep=p.k_rad_per_m<=3 & power>1e-28;
            n=n+1; expected(n)=struct('x',p.k_rad_per_m(keep),'y',power(keep)/(2*pi/160));
        end
    end
end
n2=local_match_fig(fullfile(summary,'surfaces_spectra.fig'),expected);
receipt=struct('status','ORIGINAL_FIG_DATA_PASS','failed_checks',0,'curves',n1+n2, ...
    'statistics_curves',n1,'surface_spectrum_curves',n2, ...
    'original_pixels_rerendered',false,'policy','Compare source data, not Linux/Windows pixels');
pm_sensitivity_io('write',fullfile(output_dir,'figure_audit.json'),receipt);
end

function count=local_match_fig(file,expected)
f=openfig(file,'new','invisible'); guard=onCleanup(@()close(f));
lines=findall(f,'Type','line'); assert(numel(lines)==numel(expected)); used=false(numel(expected),1);
for j=1:numel(lines)
    x=lines(j).XData(:); y=lines(j).YData(:); found=false;
    for k=1:numel(expected)
        if ~used(k) && numel(x)==numel(expected(k).x) && ...
                max(abs(x-expected(k).x))<1e-11 && max(abs(y-expected(k).y))<1e-11
            used(k)=true; found=true; break
        end
    end
    assert(found,'Original FIG curve does not match saved CSV');
end
assert(all(used)); count=numel(lines);
end
