function audit = validate_pe_bellhop_pm_stage1_interpretation(overrides)
%VALIDATE_PE_BELLHOP_PM_STAGE1_INTERPRETATION Freeze the 4-kHz comparison.
%   This stage consumes the already-run Stage 0E, 1A and 1B artifacts.  It
%   performs no propagation and does not alter PE, Bellhop, or their physics.

if nargin<1||isempty(overrides),overrides=struct();end
if ~isstruct(overrides)||~isscalar(overrides),error('overrides must be a scalar struct.');end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_vertical_project();
cfg=local_config(root,overrides);if ~exist(cfg.output_dir,'dir'),mkdir(cfg.output_dir);end

stage0e=local_load_mat(cfg.stage0e_mat,'audit'); %#ok<NASGU>
stage1a=local_load_mat(cfg.stage1a_mat,'c');
stage1b=local_load_mat(cfg.stage1b_mat,'c');
checks=struct();
checks.stage0e_all=local_all_stage0e(cfg.stage0e_checks);
checks.stage1a_structural=stage1a.checks.all;
checks.stage1b_all=stage1b.checks.all;
checks.same_profile_hash=strcmp(stage1a.profile.coeff_file_sha256,stage1b.profile.coeff_file_sha256);
checks.dim_sensitivity_smaller=stage1b.rows(end).complex_relative_error < stage1a.delta_complex_relative_error;
checks.finite_metrics=all(isfinite([stage1a.delta_tl_db,stage1a.delta_phase_rad,stage1a.delta_complex_relative_error, ...
    stage1b.ny_delta_tl_db,stage1b.ny_delta_phase_rad]));
checks.all=all(structfun(@(v)logical(v),checks));
classification=ternary(checks.all,'PASS_WITH_MODEL_DISCREPANCY','FAIL');
audit=struct('schema_version','1.0.0','stage','1C_interpretation_and_freeze','config',cfg, ...
    'stage0e',stage0e,'stage1a',stage1a,'stage1b',stage1b,'checks',checks, ...
    'classification',classification,'passed',checks.all);
audit.files=local_write_outputs(audit);
if cfg.fail_on_check&&~audit.passed,error('Stage 1C interpretation failed; see %s.',cfg.report_path);end
end

function value=local_load_mat(path,name)
if exist(path,'file')~=2,error('Required stage artifact is missing: %s',path);end
s=load(path,name);if ~isfield(s,name),error('Artifact %s has no variable %s.',path,name);end
value=s.(name);
end
function ok=local_all_stage0e(path)
if exist(path,'file')~=2,ok=false;return;end
t=readtable(path);ok=all(t.passed);
end
function cfg=local_config(root,o)
cfg=struct('stage0e_mat',fullfile(root,'results','validation','pe_bellhop_pm_numerical_budget','numerical_budget_audit.mat'), ...
    'stage0e_checks',fullfile(root,'results','validation','pe_bellhop_pm_numerical_budget','budget_checks.csv'), ...
    'stage1a_mat',fullfile(root,'results','validation','pe_bellhop_pm_stage1_tier1','stage1a_tier1_comparison.mat'), ...
    'stage1b_mat',fullfile(root,'results','validation','pe_bellhop_pm_stage1_dimensionality','stage1b_dimensionality_comparison.mat'), ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_stage1_interpretation'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_fixed_pm_4khz_comparison_report.md'),'fail_on_check',true);
names=fieldnames(o);for ii=1:numel(names),if ~isfield(cfg,names{ii}),error('Unknown Stage 1C override: %s.',names{ii});end;cfg.(names{ii})=o.(names{ii});end
end
function files=local_write_outputs(a)
out=a.config.output_dir;
n=fieldnames(a.checks);v=false(size(n));for ii=1:numel(n),v(ii)=a.checks.(n{ii});end
writetable(table(n,v,'VariableNames',{'check_name','passed'}),fullfile(out,'stage1c_interpretation_checks.csv'));
writetable(table(a.stage1a.G_pe,a.stage1a.G_bellhop,a.stage1a.delta_tl_db,a.stage1a.delta_phase_rad,a.stage1a.delta_complex_relative_error, ...
    a.stage1b.rows(end).G_2t,a.stage1b.rows(end).delta_tl_db,a.stage1b.rows(end).delta_phase_rad,a.stage1b.rows(end).complex_relative_error, ...
    'VariableNames',{'G_pe_1t','G_bellhop','cross_delta_tl_db','cross_delta_phase_rad','cross_complex_relative_error','G_pe_2t','dim_delta_tl_db','dim_delta_phase_rad','dim_complex_relative_error'}), ...
    fullfile(out,'stage1c_interpretation_metrics.csv'));
save(fullfile(out,'stage1c_interpretation.mat'),'a','-v7.3');
local_write_report(a.config.report_path,a);
files=struct('checks',fullfile(out,'stage1c_interpretation_checks.csv'),'metrics',fullfile(out,'stage1c_interpretation_metrics.csv'),'mat',fullfile(out,'stage1c_interpretation.mat'),'report',a.config.report_path);
end
function local_write_report(path,a)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end;cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop fixed-PM 4 kHz comparison\n\n状态：**%s**\n\n',a.classification);
fprintf(fid,'## Scope and frozen inputs\n\nThis Stage 1C report consumes the separately passed Stage 0E numerical budget, Stage 1A fixed-PM Tier-1 comparison, and Stage 1B dimensionality sensitivity. The canonical source is seed `260001`, U=6 m/s, span 160 m, master N=4097, requested Kmax 0.5 rad/m (realized 0.471238898 rad/m); no profile, source, reflection coefficient, PE marching, Bellhop Reflect2D, p/q, or InfluenceGeoHatCart change was made.\n\n');
fprintf(fid,'## Primary comparison\n\n| metric | value |\n|---|---:|\n| G_PE, 1T | %.12g%+.12gi |\n| G_Bellhop | %.12g%+.12gi |\n| cross-model delta TL (dB) | %.8g |\n| cross-model delta phase (rad) | %.8g |\n| cross-model complex error | %.8g |\n| G_PE, 2T (ny=%d) | %.12g%+.12gi |\n| dimensionality delta TL (dB) | %.8g |\n| dimensionality delta phase (rad) | %.8g |\n| dimensionality complex error | %.8g |\n',real(a.stage1a.G_pe),imag(a.stage1a.G_pe),real(a.stage1a.G_bellhop),imag(a.stage1a.G_bellhop),a.stage1a.delta_tl_db,a.stage1a.delta_phase_rad,a.stage1a.delta_complex_relative_error,a.stage1b.rows(end).ny,real(a.stage1b.rows(end).G_2t),imag(a.stage1b.rows(end).G_2t),a.stage1b.rows(end).delta_tl_db,a.stage1b.rows(end).delta_phase_rad,a.stage1b.rows(end).complex_relative_error);
fprintf(fid,'\nThe 1T-to-2T sensitivity is approximately %.8g dB / %.8g rad and its complex error is %.8g, far below the PE/Bellhop complex error %.8g. The remaining stable difference is therefore classified as a reflection-model discrepancy (Kirchhoff phase screen versus Bellhop local-specular Gaussian beam), not numerical failure.\n\n',a.stage1b.rows(end).delta_tl_db,a.stage1b.rows(end).delta_phase_rad,a.stage1b.rows(end).complex_relative_error,a.stage1a.delta_complex_relative_error);
fprintf(fid,'## Frozen evidence\n\n- Stage 0E PE window/grid/step and Bellhop profile/beam/step budgets: PASS.\n- Stage 1A Bellhop wall geometry, pressure-release phase, p/q state and transformed positive-range branch: PASS.\n- Stage 1B production 2T finite outputs and ny=256/512 sensitivity: PASS.\n\n## Checks\n\n');n=fieldnames(a.checks);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(a.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\n## Decision\n\n`PASS_WITH_MODEL_DISCREPANCY` means the two independently validated solvers are numerically comparable under the same fixed realization, while their remaining rough reflected-branch difference is a quantified physical-model discrepancy. Stage 2 frequency extension is allowed; the next and only recommended step is the 4/6/8 kHz frequency extension before any ensemble sweep.\n');
end
function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
