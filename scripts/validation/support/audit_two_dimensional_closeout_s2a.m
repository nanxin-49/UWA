function audit = audit_two_dimensional_closeout_s2a(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2A Independent CSV/math/integral/redraw audit.
arguments
    run_dir (1,:) char
end
data_dir=fullfile(run_dir,'data');
audit_dir=fullfile(run_dir,'audit');
assert(isfolder(data_dir) && ~isfolder(audit_dir));
t=load(fullfile(data_dir,'s2a_data.mat'),'fields','metrics','checks', ...
    'diagnostics','increment_decomposition');
f=t.fields;
source=readtable(fullfile(data_dir,'source_manifest.csv'),'TextType','string');
for ii=1:height(source)
    assert(isfile(source.path(ii)));
    assert(strcmpi(local_sha256(source.path(ii)),source.sha256(ii)));
end
line=readtable(fullfile(data_dir,'line_fields.csv'),'TextType','string');
assert(height(line)==4*3*411);
names=f.case_names; methods={'PE','BIE','FullFK'};
source_max=0; metrics_max=0; g_max=0; delta_max=0;
for j=1:4
    for m=1:3
        pick=line.case_name==string(names{j}) & line.method==string(methods{m});
        q=line(pick,:); assert(height(q)==411 && isequal(q.line_index,(1:411).'));
        u=complex(q.u_real,q.u_imag);
        source_max=max(source_max,max(abs(u-f.(methods{m}){j})));
        assert(max(abs(q.x_m-f.x_m))<1e-12);
    end
end
metric=readtable(fullfile(data_dir,'line_metrics.csv'),'TextType','string');
for row=1:height(metric)
    j=find(string(names)==metric.case_name(row));
    a=f.(char(metric.candidate(row))){j};
    b=f.(char(metric.reference(row))){j};
    w=f.weights(:); w=w/sum(w);
    e2=sum(w.*abs(a-b).^2)/sum(w.*abs(b).^2);
    mr2=sum(w.*(abs(a)-abs(b)).^2)/sum(w.*abs(b).^2);
    phase=angle(a.*conj(b));
    vals=[sqrt(e2),sqrt(mr2),sqrt(sum(w.*phase.^2)), ...
        sqrt(sum(w.*(20*log10(abs(a)./abs(b))).^2))];
    saved=[metric.complex_l2(row),metric.magnitude_l2(row), ...
        metric.phase_rms_rad(row),metric.db_rms(row)];
    metrics_max=max(metrics_max,max(abs(vals-saved)));
end
d=readtable(fullfile(data_dir,'G_delta_u.csv'),'TextType','string');
assert(height(d)==3*3*411);
for j=2:4
    for m=1:3
        pick=d.case_name==string(names{j}) & d.method==string(methods{m});
        q=d(pick,:); assert(height(q)==411);
        a=f.(methods{m}){j}; b=f.(methods{m}){1};
        g_max=max(g_max,max(abs(complex(q.G_real,q.G_imag)-a./b)));
        delta_max=max(delta_max,max(abs(complex(q.delta_u_real,q.delta_u_imag)-(a-b))));
    end
end
inc=readtable(fullfile(data_dir,'increment_decomposition.csv'),'TextType','string');
increment_max=0;
for row=1:height(inc)
    j=find(string(names)==inc.case_name(row));
    ca=char(inc.candidate(row)); rb=char(inc.reference(row));
    ef=f.(ca){1}-f.(rb){1};
    ed=(f.(ca){j}-f.(ca){1})-(f.(rb){j}-f.(rb){1});
    er=f.(ca){j}-f.(rb){j};
    den=sum(f.weights.*abs(f.(rb){1}).^2);
    vals=[sum(f.weights.*abs(ef).^2),sum(f.weights.*abs(ed).^2), ...
        2*real(sum(f.weights.*ef.*conj(ed))), ...
        sum(f.weights.*abs(er).^2)]/den;
    saved=[inc.flat_squared_term(row),inc.delta_squared_term(row), ...
        inc.cross_term(row),inc.rough_squared_term(row)];
    increment_max=max(increment_max,max(abs(vals-saved)));
end
% Reintegrate three medium receivers directly from retained quadrature data.
fk=f.FK_components{3}; k=2*pi*4000/1500;
probe=[1 206 411]; integral_max=0; component_max=0;
for ii=probe
    dx=f.x_m(ii)-fk.surface_x_m(:);
    dz=3-fk.surface_z_m(:); rr=hypot(dx,dz);
    G=1i/4*besselh(0,1,k*rr);
    dG=1i*k/4*besselh(1,1,k*rr).* ...
        (dx.*fk.surface_normal(:,1)+dz.*fk.surface_normal(:,2))./rr;
    t1=-sum(fk.quadrature_weights_m(:).*G.*fk.surface_incident_dn(:));
    t2=-sum(fk.quadrature_weights_m(:).*dG.*fk.surface_incident(:));
    integral_max=max(integral_max,abs(t1+t2-fk.receiver_field(ii)));
    component_max=max(component_max,abs(t1-fk.receiver_term_incident_dn(ii)));
    component_max=max(component_max,abs(t2-fk.receiver_term_green_dn(ii)));
end
assert(source_max<1e-12 && metrics_max<1e-12 && g_max<1e-12 && ...
    delta_max<1e-12 && increment_max<1e-12 && ...
    integral_max<1e-11 && component_max<1e-11);
mkdir(audit_dir);
redraw=fullfile(audit_dir,'redraw');
render_two_dimensional_closeout_s2a(fullfile(data_dir,'plot_data.mat'),redraw);
names={'F02_flat_control','F04_rough_native_lines','F05_native_residuals'};
pixel_equal=false(3,1); fig_reopened=false(3,1);
for i=1:3
    a=imread(fullfile(run_dir,'figures',[names{i} '.png']));
    b=imread(fullfile(redraw,[names{i} '.png']));
    pixel_equal(i)=isequal(a,b);
    h=openfig(fullfile(run_dir,'figures',[names{i} '.fig']),'invisible');
    fig_reopened(i)=~isempty(findall(h,'Type','axes'));
    close(h);
end
assert(all(pixel_equal) && all(fig_reopened));
audit=struct('pass',true,'source_hashes_checked',height(source), ...
    'line_field_rows_checked',height(line),'metric_rows_checked',height(metric), ...
    'diagnostic_rows_checked',height(d), ...
    'increment_rows_checked',height(inc), ...
    'source_complex_max_abs',source_max,'metric_max_abs',metrics_max, ...
    'G_max_abs',g_max,'delta_u_max_abs',delta_max, ...
    'increment_decomposition_max_abs',increment_max, ...
    'medium_FullFK_probe_integral_max_abs',integral_max, ...
    'medium_FullFK_component_max_abs',component_max, ...
    'png_pixel_identical',all(pixel_equal), ...
    'editable_fig_reopen_pass',all(fig_reopened), ...
    'plane_control_max_error',t.checks.plane_control_max_error);
fid=fopen(fullfile(audit_dir,'independent_audit.json'),'w','n','UTF-8');
assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
end

function h=local_sha256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
while ~feof(fid), md.update(fread(fid,1024*1024,'*uint8')); end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
