function receipt = run_d0_incident_depth(package_root,mode)
%RUN_D0_INCIDENT_DEPTH Default has no field calls; formal mode remote only.
arguments
    package_root (1,:) char
    mode (1,:) char {mustBeMember(mode,{'preflight','execute'})} = 'preflight'
end
pm_sensitivity_io('absolute',package_root);
m=jsondecode(fileread(fullfile(package_root,'manifest.json')));
pm_sensitivity_io('verify',package_root,m.files);
binding=pm_sensitivity_io('sha',fullfile(package_root,'manifest.json'));
cfg=jsondecode(fileread(fullfile(package_root,'config.json')));
assert(strcmp(cfg.kind,'D0_INCIDENT_ONLY') && isequal(cfg.depths_m(:),[100;300;600;1000;1490]));
assert(cfg.source.total_width_m==.1 && cfg.source.frequency_hz==4000 && cfg.dx_m==.05);
assert(cfg.source.v0_mps==.001 && cfg.source.c_mps==1500 && cfg.source.rho_kgpm3==1025);
assert(cfg.RAM_reserve_fraction==.30 && cfg.workers==1 && cfg.threads==1);
assert(~cfg.reflection_enabled && ~cfg.bottom_enabled && cfg.water_depth_context_m==1500);
cp=jsondecode(fileread(fullfile(package_root,'inputs','k0','checkpoint.json')));
assert(strcmp(cp.status,'AUDITED_PASS') && cp.case.seed==260001 && cp.case.N==12 && cp.case.U_PM_mps==6);
pm_sensitivity_io('verify',fullfile(package_root,'inputs','k0'),cp.artifacts_sha256);
receipt=struct('status','PREFLIGHT_PASS','depths_m',cfg.depths_m,'field_calls',0, ...
    'manifest_sha256',binding,'failed_checks',0,'representation','absolute_physical', ...
    'phase_contract_version','2.0.0');
if strcmp(mode,'preflight'), disp(jsonencode(receipt)); return; end
assert(isunix && strcmp(version('-release'),'2025b') && maxNumCompThreads==1, ...
    'Formal D0 requires Ubuntu R2025b -singleCompThread');
assert(strcmp(getenv('D0_SUPERVISED_MANIFEST_SHA256'),binding),'Use the owned remote supervisor');
out=fullfile(package_root,'output'); assert(~isfolder(out),'Existing output: inspect checkpoints; no automatic repeat/resume');
mkdir(out); diary(fullfile(out,'matlab_diary.log')); diary_guard=onCleanup(@()diary('off'));
t0=tic; cpu0=cputime; steps={};
try
    pm_sensitivity_io('write',fullfile(out,'physical_config.json'),cfg);
    k0=load(fullfile(package_root,'inputs','k0','surface.mat'),'data');
    T=k0.data.samples; ii=abs(T.x_m)<=50; T=T(ii,:);
    assert(min(T.x_m)>=-50 && max(T.x_m)<=50 && ~isempty(T));
    for depth=cfg.depths_m(:).'
        tc=tic; id=sprintf('depth_%04dm',depth); folder=fullfile(out,id); mkdir(folder);
        pm_sensitivity_io('write',fullfile(folder,'started.json'),struct('depth_m',depth,'manifest_sha256',binding));
        source=cfg.source; source.z_source_m=depth;
        x=(-220:.05:220)';
        [p,px,pz,meta]=baffled_strip_piston_field_vertical(x,0,source,64);
        data=d0_incident_diagnostics(x,p,px,pz,source,cfg); data.source_meta=meta;
        save(fullfile(folder,'field.mat'),'data','cfg','-v7');
        writetable(data.points,fullfile(folder,'points.csv'));
        writetable(data.plane_wave,fullfile(folder,'plane_wave.csv'));
        pm_sensitivity_io('write',fullfile(folder,'metrics.json'),data.metrics);
        for j=1:numel(data.window_spectra)
            w=data.window_spectra{j};
            Q=table(w.kx_rad_per_m,real(w.spectrum_f),imag(w.spectrum_f),abs(w.spectrum_f), ...
                w.phase_rad,w.phase_valid,w.propagating_band,w.outside_band,'VariableNames', ...
                {'kx_rad_per_m','real_Pa_m','imag_Pa_m','abs_Pa_m','phase_rad','phase_valid','propagating_band','outside_band_leakage'});
            writetable(Q,fullfile(folder,sprintf('window_spectrum_%d.csv',j)));
        end
        q=data.continuous_spectrum;
        writetable(table(q.prop_kx_rad_per_m,q.prop_qz_rad_per_m,real(q.prop_density_f),imag(q.prop_density_f), ...
            q.prop_phase_unwrapped_rad,'VariableNames',{'kx_rad_per_m','qz_rad_per_m','real_Pa_m','imag_Pa_m','phase_rad'}),fullfile(folder,'continuous_propagating.csv'));
        writetable(table(q.ev_kx_positive_rad_per_m,q.ev_alpha_per_m,real(q.ev_density_f),imag(q.ev_density_f), ...
            'VariableNames',{'positive_kx_rad_per_m','alpha_per_m','real_Pa_m','imag_Pa_m'}),fullfile(folder,'continuous_evanescent.csv'));
        if ismember(depth,cfg.fine_depths_m)
            fine.x_m=(-220:.025:220)';
            [fine.p_f,fine.px_f,fine.pz_f]=baffled_strip_piston_field_vertical(fine.x_m,0,source,64);
            save(fullfile(folder,'fine.mat'),'fine','-v7');
        end
        rough=struct('x_m',T.x_m,'z_m',-T.embedded_elevation_m,'depth_slope',-T.embedded_slope);
        [rough.p_f,rough.px_f,rough.pz_f]=baffled_strip_piston_field_vertical(rough.x_m,rough.z_m,source,64);
        rough.dn_water_f=(-rough.depth_slope.*rough.px_f+rough.pz_f)./hypot(1,rough.depth_slope);
        rough.normal_convention='into water n=(-dz_s/dx,+1)/hypot(1,dz_s/dx); outward=-n';
        rough.points=table(rough.x_m,rough.z_m,rough.depth_slope,real(rough.p_f),imag(rough.p_f), ...
            real(rough.px_f),imag(rough.px_f),real(rough.pz_f),imag(rough.pz_f), ...
            real(rough.dn_water_f),imag(rough.dn_water_f),'VariableNames', ...
            {'x_m','z_actual_m','depth_slope','p_real_Pa','p_imag_Pa','px_real_Pa_per_m','px_imag_Pa_per_m', ...
            'pz_real_Pa_per_m','pz_imag_Pa_per_m','dn_water_real_Pa_per_m','dn_water_imag_Pa_per_m'});
        save(fullfile(folder,'rough.mat'),'rough','-v7'); writetable(rough.points,fullfile(folder,'rough.csv'));
        audit=audit_d0_incident_case(folder,cfg);
        fprintf('D0 depth=%g m audited, elapsed=%.3f s\n',depth,toc(tc));
        paths=dir(fullfile(folder,'*')); paths={paths(~[paths.isdir]).name};
        checkpoint=struct('status','AUDITED_PASS','depth_m',depth,'config',cfg,'source',source, ...
            'manifest_sha256',binding,'input_versions_sha256',m.files,'code_versions_sha256',m.files, ...
            'completed_steps',{{'aperture_field','rough_incident','independent_audit','CSV'}}, ...
            'artifacts_sha256',pm_sensitivity_io('manifest',folder,paths),'failed_or_unrun_checks',{{}}, ...
            'audit_status_and_path','audit.json','audit_failed_checks',audit.failed_checks, ...
            'resource_usage',struct('elapsed_s',toc(tc),'CPU_cumulative_s',cputime-cpu0,'RAM','session/resource.tsv'), ...
            'execution_mode','remote_single_thread_serial','matlab_release',version('-release'), ...
            'next_allowed_step','next listed D0 depth only');
        pm_sensitivity_io('write',fullfile(folder,'checkpoint.json'),checkpoint);
        steps{end+1}=id; %#ok<AGROW>
        local_checkpoint('IN_PROGRESS',{});
    end
    render_d0_incident_depth(out,cfg,'D0 INCIDENT ONLY');
    steps{end+1}='saved_figures_and_summary';
    local_checkpoint('D0_COMPLETE_INCIDENT_ONLY',{});
    receipt=jsondecode(fileread(fullfile(out,'checkpoint.json'))); disp(jsonencode(receipt));
catch err
    fprintf(2,'%s\n',getReport(err,'extended'));
    local_checkpoint('FAILED_STOPPED',{getReport(err,'extended')});
    rethrow(err)
end
    function local_checkpoint(status,failures)
        files=dir(fullfile(out,'**','*')); files=files(~[files.isdir]); names={};
        for a=1:numel(files)
            path=fullfile(files(a).folder,files(a).name);
            rel=strrep(path(numel(out)+2:end),'\','/');
            if ismember(rel,{'checkpoint.json','matlab_diary.log','checkpoint.json.pending'}), continue; end
            names{end+1}=rel; %#ok<AGROW>
        end
        cpout=struct('status',status,'completed_steps',{steps},'config',cfg,'manifest_sha256',binding, ...
            'input_versions_sha256',m.files,'code_versions_sha256',m.files, ...
            'artifacts_sha256',pm_sensitivity_io('manifest',out,names),'failed_or_unrun_checks',{failures}, ...
            'elapsed_s',toc(t0),'CPU_s',cputime-cpu0,'resource_usage','sessions/*/resource.tsv', ...
            'audit_status_and_path','depth_*/audit.json; figures/figure_audit.json', ...
            'execution_mode','remote_supervised_single_thread','matlab_release',version('-release'), ...
            'next_allowed_step','Stop; user Fetch and Codex saved-data analysis. D1/reflection not authorized');
        pm_sensitivity_io('write',fullfile(out,'checkpoint.json'),cpout);
    end
end
