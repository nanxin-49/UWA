function result = run_bellhop_native_pm_wall_vertical(cfg)
%RUN_BELLHOP_NATIVE_PM_WALL_VERTICAL Run native ATI with the shared PM profile.
% The generic native C-ATI writer is profile-based; it does not assume a
% sinusoid, so this wrapper keeps the PM validation entrypoint explicit.
result=run_bellhop_native_sinusoidal_wall_vertical(cfg);
end
