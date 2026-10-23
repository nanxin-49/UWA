function audit=audit_two_dimensional_closeout_s2b(out,figure_folder,audit_folder)
% Independent cartesian CSV arithmetic, saved MAT/CSV and redraw audit.
if nargin<2, figure_folder='figures'; end
if nargin<3, audit_folder='audit'; end
dest=fullfile(out,audit_folder); assert(~isfolder(dest)); mkdir(dest);
t=load(fullfile(out,'data/s2b_data.mat')); p=load(fullfile(out,'data/p0_data.mat'));
F=readtable(fullfile(out,'data/fields.csv'),'TextType','string');
P=readtable(fullfile(out,'data/point_metrics.csv'),'TextType','string');
M=readtable(fullfile(out,'data/array_metrics.csv'),'TextType','string');
assert(height(F)==108 && height(P)==108 && height(M)==12);
methods={'PE','BIE','FullFK'}; cases=p.fields.case_names; maxfield=0; maxmetric=0; maxpoint=0;
for c=1:4
    for m=1:3
        r=F(F.case_name==string(cases{c})&F.method==string(methods{m}),:); assert(height(r)==9);
        assert(isequal(r.receiver_id,t.grid.receiver_id) && isequal(r.x_m,t.grid.x_m) && isequal(r.z_m,t.grid.z_m));
        maxfield=max(maxfield,max(abs(complex(r.u_native_real,r.u_native_imag)-t.U(:,c,m))));
        assert(max(abs(r.u_mapped_real-r.u_native_real))==0 && max(abs(r.u_mapped_imag-r.u_native_imag))==0);
    end
end
for j=1:height(M)
    a=F(F.case_name==M.case_name(j)&F.method==M.candidate(j),:);
    b=F(F.case_name==M.case_name(j)&F.method==M.reference(j),:);
    ar=a.u_native_real; ai=a.u_native_imag; br=b.u_native_real; bi=b.u_native_imag;
    am=hypot(ar,ai); bm=hypot(br,bi); phase=atan2(ai.*br-ar.*bi,ar.*br+ai.*bi);
    sq=(ar-br).^2+(ai-bi).^2; den=sum(br.^2+bi.^2);
    mag=sum((am-bm).^2)/den; ph=sum(2*(am.*bm-(ar.*br+ai.*bi)))/den;
    vals=[sqrt(sum(sq)/den),sqrt(mag),sqrt(mean(phase.^2)),sqrt(mean((20*log10(am./bm)).^2)),mag,ph];
    expected=[M.complex_l2(j),M.magnitude_l2(j),M.phase_rms_rad(j),M.db_rms(j),M.magnitude_squared_term(j),M.phase_squared_term(j)];
    maxmetric=max([maxmetric,max(abs(vals-expected)),abs(sum(sq)/den-mag-ph)]);
    rr=P(P.case_name==M.case_name(j)&P.candidate==M.candidate(j)&P.reference==M.reference(j),:);
    maxpoint=max([maxpoint,max(abs(rr.relative_error-sqrt(sq)./bm)),max(abs(rr.phase_signed_rad-phase))]);
end
% Re-audit all 24 line metrics independently from original CSV, not helper.
L=readtable(fullfile(out,'data/line_fields_source.csv'),'TextType','string');
LM=readtable(fullfile(out,'data/line_metrics_both_weights.csv'),'TextType','string'); maxline=0;
for j=1:height(LM)
    a=sortrows(L(L.case_name==LM.case_name(j)&L.method==LM.candidate(j),:),'line_index');
    b=sortrows(L(L.case_name==LM.case_name(j)&L.method==LM.reference(j),:),'line_index');
    if LM.weight_kind(j)=="spatial_uniform_411", w=ones(411,1)/411; else, w=a.weight/sum(a.weight); end
    ar=a.u_real; ai=a.u_imag; br=b.u_real; bi=b.u_imag;
    am=hypot(ar,ai); bm=hypot(br,bi); dp=atan2(ai.*br-ar.*bi,ar.*br+ai.*bi); den=sum(w.*bm.^2);
    mag=sum(w.*(am-bm).^2)/den; ph=sum(w.*2.*(am.*bm-ar.*br-ai.*bi))/den;
    e=sqrt(sum(w.*((ar-br).^2+(ai-bi).^2))/den);
    maxline=max([maxline,abs(e-LM.complex_l2(j)),abs(sqrt(sum(w.*dp.^2))-LM.phase_rms_rad(j)), ...
        abs(sqrt(sum(w.*(20*log10(am./bm)).^2))-LM.db_rms(j)),abs(e^2-mag-ph)]);
end
D=readtable(fullfile(out,'data/diagnostics.csv'),'TextType','string'); maxdiag=0;
for j=1:height(D)
    r=F(F.case_name==D.case_name(j)&F.method==D.method(j)&F.receiver_id==D.receiver_id(j),:);
    b=F(F.case_name=="flat"&F.method==D.method(j)&F.receiver_id==D.receiver_id(j),:);
    u=complex(r.u_native_real,r.u_native_imag); flat=complex(b.u_native_real,b.u_native_imag);
    maxdiag=max([maxdiag,abs(complex(D.G_real(j),D.G_imag(j))*flat-u),abs(flat+complex(D.delta_real(j),D.delta_imag(j))-u)]);
end
% Independent cross-term recomputation on original complex CSV.
C=readtable(fullfile(out,'data/increment.csv'),'TextType','string'); maxcross=0;
for j=1:height(C)
    ar=fieldcsv(F,C.case_name(j),C.candidate(j)); br=fieldcsv(F,C.case_name(j),C.reference(j));
    af=fieldcsv(F,"flat",C.candidate(j)); bf=fieldcsv(F,"flat",C.reference(j));
    ef=af-bf; ed=(ar-af)-(br-bf); den=sum(abs(bf).^2);
    vals=[sum(abs(ef).^2),sum(abs(ed).^2),2*real(sum(ef.*conj(ed))),sum(abs(ar-br).^2)]/den;
    maxcross=max([maxcross,max(abs(vals-[C.flat_term(j),C.delta_term(j),C.cross_term(j),C.rough_term(j)])),abs(vals(4)-sum(vals(1:3)))]);
end
% Validate raw sensitivity CSV against every corresponding MAT checkpoint.
S=readtable(fullfile(out,'data/sensitivity_fields.csv'),'TextType','string'); maxsens=0;
for job=unique(S.job_id).'
    cp=load(fullfile(out,'checkpoints',[char(job) '.mat']),'payload'); d=cp.payload;
    for variant=unique(S.variant(S.job_id==job)).'
        rr=S(S.job_id==job&S.variant==variant,:); u=complex(rr.u_real,rr.u_imag);
        if variant=="primary", v=d.u;
        elseif variant=="half_step_line411", v=d.step_half_line;
        elseif variant=="N1968_same_source", v=d.grid_double_u;
        elseif startsWith(variant,'L')
            vals=sscanf(char(variant),'L%d_N%d'); ix=cellfun(@(r) r.surface_x_m(end)==vals(1)&&numel(r.surface_x_m)==vals(2),d.runs); v=d.runs{ix}.receiver_field;
        else
            token=regexp(char(variant),'^(.*)_([0-9]+)$','tokens','once'); vv=d.level_fields.(token{1}); v=vv(:,str2double(token{2}));
        end
        maxsens=max(maxsens,max(abs(u-v)));
    end
end
% Independently recompute saved BIE spatial/window/self-quadrature changes.
R=readtable(fullfile(out,'data/reference.csv'),'TextType','string'); maxref=0;
for c=2:4
    id="bie"+string(c-1); ss=S(S.job_id==id,:); last=3+(c==4);
    fine=variantcsv(ss,"spatial_runs_"+string(last)); coarse=variantcsv(ss,"spatial_runs_"+string(last-1));
    w1=variantcsv(ss,"window_runs_1"); w2=variantcsv(ss,"window_runs_2"); w3=variantcsv(ss,"window_runs_3");
    rs=R(R.case_name==string(cases{c}),:); spatial=abs(fine-coarse); window=max([abs(w2-w1),abs(w3-w2)],[],2);
    maxref=max([maxref,max(abs(rs.spatial_abs_change-spatial)),max(abs(rs.window_abs_change-window))]);
    if c==2
        qa=variantcsv(ss,"quadrature_runs_1"); qb=variantcsv(ss,"quadrature_runs_2"); quad=abs(qa-qb);
        maxref=max(maxref,max(abs(rs.quadrature_abs_change-quad)));
    else, assert(all(isnan(rs.quadrature_abs_change))); quad=zeros(size(spatial));
    end
    maxref=max(maxref,max(abs(rs.known_component_envelope-spatial-window-quad)));
end
hashes=readtable(fullfile(out,'data/input_hashes.csv'),'TextType','string');
for j=1:height(hashes), assert(strcmpi(closeout_s2b_hash(hashes.path(j)),hashes.sha256(j))); end
assert(max([maxfield,maxmetric,maxpoint,maxline,maxdiag,maxcross,maxsens,maxref])<1e-12);
render_two_dimensional_closeout_s2b(fullfile(out,'data/plot_data.mat'),fullfile(dest,'redraw'));
png=dir(fullfile(out,figure_folder,'*.png')); assert(numel(png)==7); pixels=zeros(7,1); reopened=false(7,1);
for j=1:7
    a=imread(fullfile(png(j).folder,png(j).name)); b=imread(fullfile(dest,'redraw',png(j).name));
    assert(isequal(size(a),size(b))); pixels(j)=max(abs(double(a(:))-double(b(:))));
    [~,stem]=fileparts(png(j).name); fig=openfig(fullfile(png(j).folder,[stem '.fig']),'new','invisible');
    ax=findall(fig,'Type','axes'); assert(~isempty(ax)); reopened(j)=true; close(fig);
end
audit=struct('pass',true,'field_rows',height(F),'point_rows',height(P),'array_metric_rows',height(M),'line_metric_rows',height(LM), ...
    'sensitivity_field_rows',height(S),'source_hashes_checked',height(hashes), ...
    'max_field_csv_mat',maxfield,'max_metric_independent',maxmetric,'max_point_independent',maxpoint, ...
    'max_line_independent',maxline,'max_G_delta_identity',maxdiag,'max_cross_identity',maxcross, ...
    'max_sensitivity_csv_mat',maxsens,'max_BIE_sensitivity_csv_recompute',maxref,'png_pixel_max_difference',max(pixels),'fig_reopen_all',all(reopened), ...
    'physical_computations_during_audit',0,'physical_reference_qualified',false);
assert(all(pixels==0)); closeout_s2b_json(fullfile(dest,'independent_audit.json'),audit); disp(audit);
end
function u=fieldcsv(F,c,m)
r=F(F.case_name==c&F.method==m,:); u=complex(r.u_native_real,r.u_native_imag);
end
function u=variantcsv(S,variant)
r=sortrows(S(S.variant==variant,:),'probe_index'); u=complex(r.u_real,r.u_imag);
end
