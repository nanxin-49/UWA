function [status,reason]=classify_deepsea3d_amplitude(reference_ready,fk_grid,fk_gauss,metrics)
%CLASSIFY_DEEPSEA3D_AMPLITUDE Preserve reference/numerical/physical separation.
if ~reference_ready
    status="INCONCLUSIVE"; reason="Independent reference unresolved";
elseif ~isfinite(fk_grid) || ~isfinite(fk_gauss) || fk_grid>=1e-3 || fk_gauss>=1e-3
    status="INCONCLUSIVE"; reason="FullFK numerical guards unresolved";
elseif ~isfinite(metrics.reflected_complex) || ~isfinite(metrics.increment_complex)
    status="INCONCLUSIVE"; reason="Undefined physical comparison";
elseif metrics.reflected_complex>=.01 || metrics.increment_complex>=.01
    status="CONDITIONAL_PHYSICAL_FAIL"; reason="Frozen one-percent gate exceeded";
else
    status="CONDITIONAL_G1_PASS"; reason="";
end
end
