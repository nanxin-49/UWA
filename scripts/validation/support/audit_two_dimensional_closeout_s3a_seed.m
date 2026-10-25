function audit = audit_two_dimensional_closeout_s3a_seed(root,run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S3A_SEED Independent saved-field/CSV audit.
arguments
    root (1,:) char
    run_dir (1,:) char
end
assert(isfolder(root) && startsWith(run_dir,root));
dest=fullfile(run_dir,'audit'); assert(isfolder(dest));
d=load(fullfile(run_dir,'input','reference_data.mat'),'x','z','ix9','U','flat_exact');
ids={'pe','main96','fullfk'}; n=417; x=d.x(:); z=d.z(:); nine=d.ix9(:);
raw=readtable(fullfile(run_dir,'analysis_v2','raw_fields.csv'),TextType='string');
point=readtable(fullfile(run_dir,'analysis_v2','point_diagnostics.csv'),TextType='string');
metric=readtable(fullfile(run_dir,'analysis_v2','model_metrics.csv'),TextType='string');
phys=readtable(fullfile(run_dir,'analysis_v2','BIE_physical_response.csv'),TextType='string');
flatpair=readtable(fullfile(run_dir,'analysis_v2','flat_controls.csv'),TextType='string');
assert(height(raw)==3*n && height(point)==2*n && height(metric)==4);
u=complex(zeros(n,3)); f=squeeze(d.U(:,1,:));
for m=1:3
    j=load(fullfile(run_dir,'jobs',[ids{m} '.mat']),'payload','receipt');
    assert(strcmp(j.receipt.status,'COMPLETE') && j.receipt.primary_points==n);
    assert(isequal(j.payload.x_m(:),x) && isequal(j.payload.z_m(:),z));
    u(:,m)=j.payload.u_native(:);
end
maxraw=0;
for m=1:3
    ii=(m-1)*n+(1:n); r=raw(ii,:);
    expected=[real(u(:,m)),imag(u(:,m)),real(f(:,m)),imag(f(:,m))];
    actual=[r.u_rough_real,r.u_rough_imag,r.u_flat_real,r.u_flat_imag];
    maxraw=max(maxraw,max(abs(expected-actual),[],'all'));
    assert(isequal(r.x_m,x) && isequal(r.z_m,z));
    assert(max(abs((r.delta_u_real+1i*r.delta_u_imag)-(u(:,m)-f(:,m))))<1e-14);
end
assert(maxraw<1e-14);
b=load(fullfile(run_dir,'jobs','main96.mat'),'payload'); b=b.payload.u_native(:);
maxmetric=0; maxpoint=0;
for m=[1 3]
    if m==1, name="PE"; j=1; else, name="FullFK"; j=2; end
    rows=(j-1)*n+(1:n); rr=point(rows,:); e=u(:,m)-b;
    maxpoint=max(maxpoint,max(abs(rr.residual_real+1i*rr.residual_imag-e)));
    assert(max(abs(rr.residual_abs-abs(e)))<1e-14);
    assert(max(abs(rr.wrapped_phase_difference_rad-angle(u(:,m).*conj(b))))<1e-14);
    assert(max(abs(abs(e).^2-(abs(u(:,m))-abs(b)).^2 ...
        -4*abs(u(:,m)).*abs(b).*sin(angle(u(:,m).*conj(b))/2).^2))<1e-12);
    for s=1:2
        if s==1, ii=nine; sid="nine_9"; else, ii=(1:411)'; sid="line_411"; end
        row=metric(metric.method==name & metric.support_id==sid,:);
        assert(height(row)==1 && row.n==numel(ii));
        recompute=norm(e(ii))/norm(b(ii));
        maxmetric=max(maxmetric,abs(row.complex_L2-recompute));
        [v,k]=max(abs(e(ii))); assert(abs(v-row.max_abs_residual)<1e-14);
        assert(ii(k)==row.worst_index && x(ii(k))==row.worst_x_m && z(ii(k))==row.worst_z_m);
        assert(abs(row.BIE_rough_denominator_sum_abs2-sum(abs(b(ii)).^2))<1e-12);
    end
end
assert(maxmetric<1e-14 && maxpoint<1e-14);
for s=1:2
    if s==1, ii=nine; sid="nine_9"; else, ii=(1:411)'; sid="line_411"; end
    row=phys(phys.support_id==sid,:); assert(height(row)==1);
    assert(abs(row.BIE_physical_response_L2-norm(b(ii)-f(ii,2))/norm(f(ii,2)))<1e-14);
    for m=[1 3]
        if m==1, name="PE"; else, name="FullFK"; end
        q=flatpair(flatpair.method==name & flatpair.support_id==sid,:);
        assert(height(q)==1 && abs(q.flat_complex_L2_to_BIE- ...
            norm(f(ii,m)-f(ii,2))/norm(f(ii,2)))<1e-14);
    end
end
levels={'self64','space14','window60','window70'};
maxdiff=zeros(n,1);
for j=1:4
    z0=load(fullfile(run_dir,'jobs',[levels{j} '.mat']),'payload');
    maxdiff=max(maxdiff,abs(b-z0.payload.u_native(:)));
end
assert(max(abs(point.BIE_max_level_difference_abs(1:n)-maxdiff))<1e-14);
assert(max(abs(point.BIE_max_level_difference_abs(n+1:end)-maxdiff))<1e-14);
pe0=load(fullfile(run_dir,'jobs','pe.mat'),'payload');
peg=load(fullfile(run_dir,'jobs','pe_grid.mat'),'payload','receipt');
assert(strcmp(peg.receipt.status,'COMPLETE'));
pe_sensitivity=abs(u(:,1)-pe0.payload.half_step_native(:)) ...
    +abs(u(:,1)-peg.payload.u_native(:));
assert(max(abs(point.PE_sensitivity_abs(1:n)-pe_sensitivity))<1e-14);
figure_names={'F01_PM_geometry_v2','F02_nine_point_errors_v2','F03_method_line_v2'};
figure_ok=false(3,1); png_bytes=zeros(3,1); pixel_differences=zeros(3,1);
for j=1:3
    png=fullfile(run_dir,'figures',[figure_names{j} '.png']);
    fig=fullfile(run_dir,'figures',[figure_names{j} '.fig']);
    info=imfinfo(png); assert(info.Width>=1000 && info.Height>=600);
    h=openfig(fig,'invisible');
    assert(~isempty(findall(h,'Type','axes')));
    redraw=fullfile(dest,[figure_names{j} '_fig_reexport_v3.png']);
    assert(~isfile(redraw)); drawnow;
    exportgraphics(h,redraw,'Resolution',180); close(h);
    original_pixels=imread(png); redraw_pixels=imread(redraw);
    assert(isequal(size(original_pixels),size(redraw_pixels)));
    pixel_differences(j)=nnz(original_pixels~=redraw_pixels);
    png_bytes(j)=dir(png).bytes; figure_ok(j)=true;
end
audit=struct('status','S3A_SEED260001_SAVED_DATA_AUDIT_PASS', ...
    'raw_rows',height(raw),'point_rows',height(point),'metric_rows',height(metric), ...
    'nine_point_count',numel(nine),'full_line_count',411, ...
    'max_raw_complex_difference',maxraw,'max_point_complex_difference',maxpoint, ...
    'max_metric_difference',maxmetric,'BIE_five_level_field_difference_max',max(maxdiff), ...
    'figure_fig_reopen',figure_ok,'png_bytes',png_bytes, ...
    'fig_reexport_pixel_differences',pixel_differences, ...
    'figure_pixel_redraw_status','FIG object reexport compared; independent plotting source not executed');
fid=fopen(fullfile(dest,'independent_saved_data_audit_v3.json'),'w');
assert(fid>0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); disp(audit);
end
