function audit=audit_two_dimensional_closeout_s2b_reference(out)
% Independent saved-CSV audit of current reference-qualified field candidate.
dest=fullfile(out,'formal_reference_v2'); ad=fullfile(dest,'audit'); assert(~isfolder(ad)); mkdir(ad);
d=load(fullfile(dest,'reference_data.mat')); p=load(fullfile(out,'data/p0_data.mat'));
core=readtable(fullfile(dest,'fields_core.csv'),'TextType','string'); line=readtable(fullfile(dest,'fields_line.csv'),'TextType','string');
lm=readtable(fullfile(dest,'line_metrics.csv'),'TextType','string'); am=readtable(fullfile(dest,'array_metrics.csv'),'TextType','string');
assert(height(core)==108&&height(line)==4932&&height(lm)==24&&height(am)==12);
maxfield=0; maxmetric=0; maxpoint=0; maxidentity=0; methods={'PE','BIE','FullFK'}; cases=p.fields.case_names;
for c=1:4
    for m=1:3
        a=raw(line,cases{c},methods{m}); b=raw(core,cases{c},methods{m});
        maxfield=max([maxfield,max(abs(a-d.U(1:411,c,m))),max(abs(b-d.U(d.ix9,c,m)))]);
        if m<3
            assert(max(abs(d.U(:,c,m)-d.U_legacy(:,c,m)))==0,'PE/BIE primary was altered');
        end
    end
end
for j=1:height(lm)
    a=raw(line,lm.case_name(j),lm.candidate(j)); b=raw(line,lm.case_name(j),lm.reference(j));
    if lm.weight_kind(j)=="spatial_uniform_411", w=ones(411,1)/411; else, w=p.fields.weights; end
    [vals,identity]=cartesian(a,b,w);
    expected=[lm.complex_l2(j),lm.magnitude_l2(j),lm.phase_rms_rad(j),lm.db_rms(j),lm.magnitude_squared_term(j),lm.phase_squared_term(j)];
    maxmetric=max(maxmetric,max(abs(vals-expected))); maxidentity=max(maxidentity,identity);
end
pt=readtable(fullfile(dest,'point_metrics.csv'),'TextType','string');
for j=1:height(am)
    a=raw(core,am.case_name(j),am.candidate(j)); b=raw(core,am.case_name(j),am.reference(j));
    [vals,identity]=cartesian(a,b,ones(9,1)/9);
    expected=[am.complex_l2(j),am.magnitude_l2(j),am.phase_rms_rad(j),am.db_rms(j),am.magnitude_squared_term(j),am.phase_squared_term(j)];
    maxmetric=max(maxmetric,max(abs(vals-expected))); maxidentity=max(maxidentity,identity);
    rr=pt(pt.case_name==am.case_name(j)&pt.candidate==am.candidate(j)&pt.reference==am.reference(j),:);
    maxpoint=max([maxpoint,max(abs(rr.relative_error-abs(a-b)./abs(b))),max(abs(rr.denominator_squared-abs(b).^2))]);
end
S=readtable(fullfile(dest,'sensitivity_fields.csv'),'TextType','string'); R=readtable(fullfile(dest,'reference.csv'),'TextType','string');
maxref=0; maxnumerical=0;
for c=1:4
    hi=variant(S,c,'FullFK_complete_GL800'); lo=variant(S,c,'FullFK_complete_GL400');
    tail=variant(S,c,'FullFK_analytic_exterior_trapezoid4097');
    maxnumerical=max([maxnumerical,abs(norm(hi-lo)/norm(hi)-d.numerical.FullFK_GL400_to800_relative(c)), ...
        abs(max(abs(hi-tail))-d.numerical.FullFK_independent_tail_max_abs(c))]);
    if c>1
        u64=variant(S,c,'BIE_q64_current'); u96=variant(S,c,'BIE_q96'); sp=variant(S,c,'BIE_spatial_previous');
        w1=variant(S,c,'BIE_window60'); w2=variant(S,c,'BIE_window70'); w3=variant(S,c,'BIE_window80');
        estimate=abs(u96-sp)+max(abs([w2-w1,w3-w2]),[],2)+abs(u64-u96);
        r=R(R.case_name==string(cases{c}),:); maxref=max(maxref,max(abs(estimate-r.resolution_estimate)));
    end
    for j=1:2
        run=d.fk_runs{c,j}; assert(isinf(run.local_radius_m));
        maxidentity=max(maxidentity,max(abs(run.receiver_field-run.term_incident_dn-run.term_green_dn)));
        [eta,ep]=closeout_s2b_surface(run.surface_x_m,[0 .01 .05 .2]*(1:4==c)');
        assert(max(abs(eta-run.surface_z_m))<1e-14);
        assert(max(abs(run.surface_normal-[-ep./hypot(1,ep),1./hypot(1,ep)]),[],'all')<1e-14);
    end
end
G=readtable(fullfile(dest,'diagnostics.csv'),'TextType','string'); C=readtable(fullfile(dest,'increments.csv'),'TextType','string');
maxdiag=0; maxcross=0;
for j=1:height(G)
    a=raw_at(core,G.case_name(j),G.method(j),G.receiver_id(j)); b=raw_at(core,'flat',G.method(j),G.receiver_id(j));
    maxdiag=max([maxdiag,abs(complex(G.G_real(j),G.G_imag(j))*b-a),abs(b+complex(G.delta_real(j),G.delta_imag(j))-a)]);
end
for j=1:height(C)
    ar=raw(core,C.case_name(j),C.candidate(j)); br=raw(core,C.case_name(j),C.reference(j));
    af=raw(core,'flat',C.candidate(j)); bf=raw(core,'flat',C.reference(j));
    ef=af-bf; ed=(ar-af)-(br-bf); den=sum(abs(bf).^2);
    vv=[sum(abs(ef).^2),sum(abs(ed).^2),2*real(sum(ef.*conj(ed))),sum(abs(ar-br).^2)]/den;
    maxcross=max([maxcross,max(abs(vv-[C.flat_term(j),C.delta_term(j),C.cross_term(j),C.rough_term(j)])),abs(vv(4)-sum(vv(1:3)))]);
end
% Check the original P0 certificate and immutable source/code chain again.
assert(p.audit.pass && p.audit.physical_field_computations==0);
H=readtable(fullfile(out,'data/input_hashes.csv'),'TextType','string');
for j=1:height(H), assert(strcmpi(closeout_s2b_hash(H.path(j)),H.sha256(j))); end
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
solver_hash=closeout_s2b_hash(fullfile(root,'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m'));
for j=1:numel(d.bcost), assert(strcmp(d.bcost{j}.solver_sha256,solver_hash) && d.bcost{j}.wall_seconds<900); end
assert(all(d.density_checks.current_saved_density_relative<1e-10));
assert(max([maxfield,maxmetric,maxpoint,maxidentity,maxref,maxnumerical,maxdiag,maxcross])<1e-12);
assert(d.summary.charged_total_seconds<7200);
render_two_dimensional_closeout_s2b(fullfile(dest,'plot_data.mat'),fullfile(dest,'figures'));
render_two_dimensional_closeout_s2b(fullfile(dest,'plot_data.mat'),fullfile(ad,'redraw'));
F=dir(fullfile(dest,'figures/*.png')); assert(numel(F)==7); pixels=0;
for j=1:7
    a=imread(fullfile(F(j).folder,F(j).name)); b=imread(fullfile(ad,'redraw',F(j).name)); assert(isequal(size(a),size(b)));
    pixels=max(pixels,max(abs(double(a(:))-double(b(:)))));
    [~,name]=fileparts(F(j).name); fig=openfig(fullfile(dest,'figures',[name '.fig']),'new','invisible');
    assert(~isempty(findall(fig,'Type','axes'))); close(fig);
end
assert(pixels==0);
audit=struct('pass',true,'core_field_rows',height(core),'line_field_rows',height(line),'line_metric_rows',height(lm),'array_metric_rows',height(am), ...
    'sensitivity_field_rows',height(S),'max_csv_mat_field',maxfield,'max_independent_metric',maxmetric, ...
    'max_point_denominator',maxpoint,'max_identity',maxidentity,'max_BIE_resolution_csv',maxref, ...
    'max_FullFK_independent_numerical_csv',maxnumerical,'max_G_delta',maxdiag,'max_cross',maxcross, ...
    'png_pixel_max_difference',pixels,'fig_reopen_all',true,'original_P0_preserved',true,'input_hashes_rechecked',height(H), ...
    'new_physics_during_audit',0,'reference_interpretation','Qualified for the explicitly reported empirical numerical resolution, not a rigorous error bound or model physical PASS');
closeout_s2b_json(fullfile(ad,'independent_audit.json'),audit); disp(audit);
end
function u=raw(T,c,m)
r=T(T.case_name==string(c)&T.method==string(m),:); u=complex(r.u_real,r.u_imag);
end
function u=raw_at(T,c,m,id)
r=T(T.case_name==string(c)&T.method==string(m)&T.receiver_id==string(id),:); u=complex(r.u_real,r.u_imag);
end
function u=variant(S,c,name)
r=sortrows(S(S.case_index==c&S.variant==string(name),:),'sample_index'); u=complex(r.u_real,r.u_imag);
end
function [v,identity]=cartesian(a,b,w)
ar=real(a); ai=imag(a); br=real(b); bi=imag(b); aa=hypot(ar,ai); bb=hypot(br,bi);
sq=(ar-br).^2+(ai-bi).^2; den=sum(w.*bb.^2); phase=atan2(ai.*br-ar.*bi,ar.*br+ai.*bi);
mag=sum(w.*(aa-bb).^2)/den; ph=sum(w.*2.*(aa.*bb-ar.*br-ai.*bi))/den;
v=[sqrt(sum(w.*sq)/den),sqrt(mag),sqrt(sum(w.*phase.^2)),sqrt(sum(w.*(20*log10(aa./bb)).^2)),mag,ph];
identity=abs(sum(w.*sq)/den-mag-ph);
end
