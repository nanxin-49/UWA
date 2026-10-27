function result=audit_deepsea3d_65grid_four_figures(project_root)
%AUDIT_DEEPSEA3D_65GRID_FOUR_FIGURES Verify saved FIG objects against native CSV.
arguments
    project_root (1,1) string
end
out=fullfile(project_root,'results','validation','deepsea3d_65grid_four_20261027');
t=readtable(fullfile(out,'four_point_blind_metrics.csv'));
f1=fullfile(out,'four_point_reflection_complex_v3.fig');
f2=fullfile(out,'prediction_vs_FE_sensitivity_v3.fig');
assert(isfile(f1) && isfile(f2));
a=openfig(f1,'invisible');
objects=findall(a,'Type','line');
for j=1:height(t)
    expected=[t.R_FE_real(j)*1e5,t.R_FE_imag(j)*1e5];
    residual=[(t.Rpred_real(j)-t.R_FE_real(j))*1e9, ...
        (t.Rpred_imag(j)-t.R_FE_imag(j))*1e9];
    assert(local_has_point(objects,expected,1e-10));
    assert(local_has_point(objects,residual,1e-10));
end
close(a);
b=openfig(f2,'invisible');
objects=findall(b,'Type','line');
assert(local_has_series(objects,t.absolute_complex_gap_m_inv,1e-18));
assert(local_has_series(objects,t.FE_empirical_sensitivity_m_inv,1e-18));
assert(local_has_series(objects,5*t.FE_empirical_sensitivity_m_inv,1e-18));
assert(all(arrayfun(@(h)isempty(h.LineStyle) || strcmp(h.LineStyle,'none'),objects)));
close(b);
result=struct('status','PASS','figures_reopened',2,'four_native_points_verified',height(t), ...
    'residual_points_verified',height(t),'three_error_series_verified',true, ...
    'discrete_marker_only',true);
fid=fopen(fullfile(out,'figure_object_audit.json'),'w');assert(fid>0);
fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
disp(result);
end

function tf=local_has_point(objects,p,tol)
tf=false;
for h=reshape(objects,1,[])
    x=h.XData;y=h.YData;
    if isscalar(x) && isscalar(y) && isfinite(x) && isfinite(y) && ...
            max(abs([x y]-p))<tol
        tf=true;return
    end
end
end

function tf=local_has_series(objects,values,tol)
tf=false;
for h=reshape(objects,1,[])
    y=h.YData;
    if numel(y)==numel(values) && max(abs(y(:)-values(:)))<tol
        tf=true;return
    end
end
end
