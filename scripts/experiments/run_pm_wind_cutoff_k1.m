function plan = run_pm_wind_cutoff_k1(varargin)
%RUN_PM_WIND_CUTOFF_K1 Plan-only entry. No flag can enable an acoustic solve.
plan=struct('status','DISABLED_PENDING_K0_REVIEW_AND_EXPLICIT_K1_AUTHORIZATION', ...
    'seed',260001,'provisional_U_N',[6 12;6 7;6 40;5 12;8 12;5 40;8 40], ...
    'new_BIE_per_surface',true,'receiver_count',411,'frequency_hz',4000,'source_width_m',.10);
if nargin>0
    error('PM:K1Disabled','K1 execution is not implemented or authorized. Review K0 and approve an explicit numerical/resource policy first.');
end
disp(jsonencode(plan));
end
