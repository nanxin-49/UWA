function audit = audit_two_dimensional_closeout_s0_s1(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S0_S1 Independent CSV/MAT and image replay.
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
assert(isfolder(run_dir));
t=load(fullfile(run_dir,'plot_data','plot_data.mat'),'plot_data'); d=t.plot_data;
line=readtable(fullfile(run_dir,'line_metrics.csv'),'TextType','string');
line_points=readtable(fullfile(run_dir,'line_point_metrics.csv'),'TextType','string');
points=readtable(fullfile(run_dir,'point_metrics.csv'),'TextType','string');
summaries=readtable(fullfile(run_dir,'point_summary.csv'),'TextType','string');
native=readtable(fullfile(run_dir,'native_field_samples.csv'),'TextType','string');
names={d.line.cases.name}; pair_names=["PE_BIE","BH_BIE","PE_BH"];
csv_paths={"weak_low_K.csv","region_II_low_K.csv","strong_height_low_K.csv"};
line_max=0;
g_point_max=0;
for ci=1:3
    p=fullfile(root,'results','visualization','pe_bellhop_bie_comparison', ...
        '20261018_194331_523','plot_data',csv_paths{ci});
    source=readtable(p);
    ix=logical(source.M99_mask);
    w=source.incident_weight(ix); w=w/sum(w);
    fields={complex(source.PE_real,source.PE_imag), ...
        complex(source.BH_real,source.BH_imag),complex(source.BIE_real,source.BIE_imag)};
    for j=1:height(points)
        q=points(j,:);
        if q.case_id~=string(names{ci}) || q.quantity_kind~="G" || ~q.valid, continue; end
        row=find(abs(source.x_m-q.x_m)<=1e-9); assert(isscalar(row) && q.z_m==3);
        pair=find(pair_names==q.pair); assert(isscalar(pair));
        l=[1 2 1]; r=[3 3 2];
        source_A=fields{l(pair)}; source_B=fields{r(pair)};
        value_A=source_A(row); value_B=source_B(row);
        g_point_max=max([g_point_max,abs(value_A-complex(q.A_real,q.A_imag)), ...
            abs(value_B-complex(q.B_real,q.B_imag))]);
    end
    left=[1 2 1]; right=[3 3 2];
    for pi=1:3
        A=fields{left(pi)}; B=fields{right(pi)}; A=A(ix); B=B(ix);
        v=[sqrt(sum(w.*abs(A-B).^2)/sum(w.*abs(B).^2)), ...
            sqrt(sum(w.*angle(A.*conj(B)).^2)), ...
            sqrt(sum(w.*(20*log10(abs(A)./abs(B))).^2))];
        row=line(line.case_id==string(names{ci})&line.pair==pair_names(pi),:);
        assert(height(row)==1 && row.n_valid==nnz(ix));
        target=[row.complex_l2,row.phase_rms_rad,row.magnitude_ratio_rms_db];
        line_max=max(line_max,max(abs(v-target)));
        assert(all(abs(v-target)<=1e-12+1e-10*abs(target)), ...
            'CSV line rebuild differs from saved MAT line metrics.');
        point_table=line_points(line_points.case_id==string(names{ci}) & ...
            line_points.pair==pair_names(pi),:);
        assert(height(point_table)==nnz(ix));
        assert(max(abs(point_table.x_m-source.x_m(ix)))<1e-9);
        assert(max(abs(point_table.absolute_difference-abs(A-B)))<1e-12);
        assert(max(abs(point_table.relative_difference-abs(A-B)./abs(B)))<1e-12);
    end
end
assert(g_point_max<1e-12,'Point G values disagree with source CSV.');
raw_csv=readtable(fullfile(root,'results','visualization','pe_bellhop_bie_xz', ...
    '20261019_115214_348','plot_data','xz_fields.csv'));
raw_max=0;
for j=1:height(native)
    q=native(j,:);
    ix=find(abs(raw_csv.x_m-q.x_m)<=1e-9 & abs(raw_csv.z_m-q.z_m)<=1e-9);
    assert(isscalar(ix));
    switch q.method
        case "PE"
            value=complex(raw_csv.PE_carrier_restored_real(ix),raw_csv.PE_carrier_restored_imag(ix));
        case "BIE"
            value=complex(raw_csv.BIE_native_real(ix),raw_csv.BIE_native_imag(ix));
        case "BH"
            value=complex(raw_csv.BH_native_real(ix),raw_csv.BH_native_imag(ix));
    end
    raw_max=max(raw_max,abs(value-complex(q.real,q.imag)));
end
assert(raw_max<1e-12,'Native raw point values disagree with source CSV.');
point_max=0;
for j=1:height(points)
    q=points(j,:);
    if ~q.valid, assert(isnan(q.complex_abs)); continue; end
    A=complex(q.A_real,q.A_imag); B=complex(q.B_real,q.B_imag);
    delta=abs(A-B);
    point_max=max(point_max,abs(delta-q.complex_abs));
    assert(abs(delta-q.complex_abs)<1e-12);
    if abs(B)>0
        assert(abs(delta/abs(B)-q.relative_abs)<1e-12);
        assert(abs(20*log10(abs(A)/abs(B))-q.magnitude_db)<1e-12);
    end
    if q.quantity_kind~="raw_native_amplitude_only"
        assert(abs(angle(A*conj(B))-q.phase_rad)<1e-12);
        assert(abs(q.complex_abs^2-q.magnitude_sq-q.phase_sq)<1e-12);
    else
        assert(isnan(q.phase_rad));
    end
end
for j=1:height(summaries)
    s=summaries(j,:);
    q=points(points.case_id==s.case_id & points.pair==s.pair & ...
        points.quantity_kind==s.quantity_kind,:);
    assert(height(q)==9 && nnz(q.valid)==s.n_valid);
    assert(nnz(isfinite(q.relative_abs))==s.n_relative);
    assert(nnz(isfinite(q.magnitude_db))==s.n_db);
    assert(nnz(isfinite(q.phase_rad))==s.n_phase);
    if s.n_valid>0
        active=logical(q.valid);
        A=complex(q.A_real(active),q.A_imag(active));
        B=complex(q.B_real(active),q.B_imag(active));
        val=norm(A-B)/norm(B);
        assert(abs(val-s.normalized_error)<=1e-12+1e-10*abs(val));
    end
end
figs=["F01","F03","F04","F05","F05_grid","F06","F07"];
same=false(numel(figs),1);
for j=1:numel(figs)
    a=imread(fullfile(run_dir,'figures',figs(j)+".png"));
    b=imread(fullfile(run_dir,'audit','redraw',figs(j)+".png"));
    same(j)=isequal(a,b);
end
assert(all(same));
audit=struct('csv_line_max_abs_difference',line_max, ...
    'g_point_source_csv_max_abs_difference',g_point_max, ...
    'raw_point_source_csv_max_abs_difference',raw_max, ...
    'point_complex_max_abs_difference',point_max, ...
    'line_pairs_checked',height(line),'line_point_rows_checked',height(line_points), ...
    'point_rows_checked',height(points), ...
    'native_field_rows_checked',height(native),'summary_rows_checked',height(summaries), ...
    'figures_pixel_identical',all(same), ...
    'new_physical_field_calculations',0,'pass',true);
fid=fopen(fullfile(run_dir,'audit','independent_audit.json'),'w','n','UTF-8');
assert(fid>=0); cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); clear cleaner
disp(audit);
end
