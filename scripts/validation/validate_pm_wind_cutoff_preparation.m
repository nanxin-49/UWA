function receipt = validate_pm_wind_cutoff_preparation(package_root, output_dir)
%VALIDATE_PM_WIND_CUTOFF_PREPARATION Local tests only; never official K0/K1.
pm_sensitivity_io('absolute',output_dir); assert(~isfolder(output_dir)); mkdir(output_dir);
diary(fullfile(output_dir,'matlab.log')); dg=onCleanup(@()diary('off')); t0=tic;
old_threads=maxNumCompThreads(1); tg=onCleanup(@()maxNumCompThreads(old_threads));
cfg=jsondecode(fileread(fullfile(package_root,'config.json'))); tests=0;
sources=jsondecode(fileread(fullfile(package_root,'provenance','sources.json')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
for j=1:numel(sources.sources)
    assert(strcmp(pm_sensitivity_io('sha',fullfile(root,sources.sources(j).path)),sources.sources(j).sha256));
end
before=rng;
recovery=zeros(3,1); max_pair_error=0; configs=0;
for seed=260001:260003
    anchor=readtable(fullfile(package_root,'inputs',sprintf('seed%d.csv',seed)));
    all=build_pm_wind_cutoff_surface(anchor,seed,6,40,256);
    recovery(seed-260000)=all.metrics.gaussian_recovery_max_error;
    assert(isequal(all.coefficients.cos_depth_coeff_m(1:12),anchor.cos_coeff_m));
    assert(isequal(all.coefficients.sin_depth_coeff_m(1:12),anchor.sin_coeff_m)); tests=tests+2;
    for U=[5 6 8]
        whole=build_pm_wind_cutoff_surface(anchor,seed,U,40,256);
        for N=[7 12 19 32 40]
            d=build_pm_wind_cutoff_surface(anchor,seed,U,N,256); configs=configs+1;
            assert(isequal(table2array(d.coefficients),table2array(whole.coefficients(1:N,:))));
            assert(isequal(d.coefficients.gaussian_cos,all.coefficients.gaussian_cos(1:N)));
            assert(isequal(d.coefficients.gaussian_sin,all.coefficients.gaussian_sin(1:N))); tests=tests+3;
            sp=raw_pm_spectrum_grid_vertical(U,256,2,160,160);
            difference=max(abs(d.coefficients.Sk_m3-sp.E1D(1,2:N+1).'));
            assert(difference<1e-12); max_pair_error=max(max_pair_error,difference); tests=tests+1;
            assert(d.metrics.Knyquist_rad_per_m>d.metrics.Kmax_rad_per_m);
            assert(abs(sum(d.spectra.raw_bin_power_m2)-d.metrics.raw_rms_height_m^2)<1e-12); tests=tests+2;
        end
    end
end
assert(isequal(before,rng)); tests=tests+1;
% Default is preflight, official execution is blocked on this Windows test.
p=run_pm_wind_cutoff_k0(package_root); assert(strcmp(p.status,'PREFLIGHT_PASS')); tests=tests+1;
plan=run_pm_wind_cutoff_k1(); assert(contains(plan.status,'DISABLED')); tests=tests+1;
local_reject(@()run_pm_wind_cutoff_k1('execute'),'K1 enabled by argument'); tests=tests+1;
% Full save/audit/resume/render integration on THREE TEST_ONLY cases.
fixture=fullfile(output_dir,'test_package'); mkdir(fixture);
copyfile(fullfile(package_root,'inputs'),fullfile(fixture,'inputs'));
copyfile(fullfile(package_root,'code'),fullfile(fixture,'code'));
cfg.kind='TEST_ONLY'; cfg.cases=struct('seed',{260001,260001,260001}, ...
    'U_PM_mps',{8,5,6},'N',{40,7,12});
pm_sensitivity_io('write',fullfile(fixture,'config.json'),cfg);
files=dir(fullfile(fixture,'**','*')); files=files(~[files.isdir]); rel=cell(1,numel(files));
for j=1:numel(files), rel{j}=strrep(fullfile(files(j).folder,files(j).name),[fixture filesep],''); end
manifest=struct('files',pm_sensitivity_io('manifest',fixture,rel));
pm_sensitivity_io('write',fullfile(fixture,'manifest.json'),manifest);
local_reject(@()run_pm_wind_cutoff_k0(fixture,'execute'),'Test package in official mode'); tests=tests+1;
r=run_pm_wind_cutoff_k0(fixture,'test-only'); assert(r.case_count==3 && r.failed_checks==0); tests=tests+1;
seal=fullfile(fixture,'output','checkpoint.json'); old_sha=pm_sensitivity_io('sha',seal);
r=run_pm_wind_cutoff_k0(fixture,'test-only'); assert(strcmp(r.action,'VERIFIED_NO_REPEAT'));
assert(strcmp(old_sha,pm_sensitivity_io('sha',seal))); tests=tests+1;
% Audit-only recovery and case-skip recovery must preserve generation bytes.
case_path=fullfile(fixture,'output','seed260001_U8_N40');
generation_sha=pm_sensitivity_io('sha',fullfile(case_path,'generation.log'));
movefile(seal,fullfile(output_dir,'aggregate_checkpoint_before_resume.json'));
movefile(fullfile(case_path,'checkpoint.json'),fullfile(output_dir,'case_checkpoint_before_reaudit.json'));
r=run_pm_wind_cutoff_k0(fixture,'test-only');
assert(r.case_count==3 && strcmp(generation_sha,pm_sensitivity_io('sha',fullfile(case_path,'generation.log')))); tests=tests+1;
% Output drift must be rejected even after aggregate completion.
victim=fullfile(fixture,'output','seed260001_U8_N40','metrics.json'); original=fileread(victim);
fid=fopen(victim,'a'); fprintf(fid,' '); fclose(fid);
local_reject(@()run_pm_wind_cutoff_k0(fixture,'test-only'),'Modified completed artifact accepted'); tests=tests+1;
fid=fopen(victim,'w'); fprintf(fid,'%s',original); fclose(fid);
% Isolated interrupted-generation fixture; never implicitly regenerate.
orphan=fullfile(output_dir,'orphan_package'); mkdir(orphan);
copyfile(fullfile(fixture,'inputs'),fullfile(orphan,'inputs'));
copyfile(fullfile(fixture,'code'),fullfile(orphan,'code'));
copyfile(fullfile(fixture,'manifest.json'),fullfile(orphan,'manifest.json'));
copyfile(fullfile(fixture,'config.json'),fullfile(orphan,'config.json'));
folder=fullfile(orphan,'output','seed260001_U8_N40'); mkdir(folder);
pm_sensitivity_io('write',fullfile(folder,'generation_started.json'),struct('interrupted_test',true));
local_reject(@()run_pm_wind_cutoff_k0(orphan,'test-only'),'Interrupted generation repeated'); tests=tests+1;
% Input drift / changed cutoff are rejected before output creation.
victim=fullfile(orphan,'inputs','seed260001.csv'); fid=fopen(victim,'a'); fprintf(fid,' '); fclose(fid);
local_reject(@()run_pm_wind_cutoff_k0(orphan),'Input drift accepted'); tests=tests+1;
for j=1:numel(sources.sources)
    assert(strcmp(pm_sensitivity_io('sha',fullfile(root,sources.sources(j).path)),sources.sources(j).sha256));
end
receipt=struct('status','LOCAL_PREPARATION_TESTS_PASS','test_assertions',tests,'coefficient_configurations_tested',configs, ...
    'saved_surface_test_cases',3,'official_K0_cases_executed',0,'acoustic_calls',0,'SSH_calls',0, ...
    'gaussian_recovery_errors',recovery,'spectrum_difference_from_production_formula',max_pair_error, ...
    'failed_checks',0,'elapsed_s',toc(t0),'matlab_release',version('-release'), ...
    'limitations',{{'Linux flock/nohup and real SSH not executed','No official 45-case K0 or K1 acoustic validation'}});
pm_sensitivity_io('write',fullfile(output_dir,'receipt.json'),receipt); disp(jsonencode(receipt));
end

function local_reject(f, why)
rejected=false;
try
    f();
catch
    rejected=true;
end
assert(rejected,why);
end
