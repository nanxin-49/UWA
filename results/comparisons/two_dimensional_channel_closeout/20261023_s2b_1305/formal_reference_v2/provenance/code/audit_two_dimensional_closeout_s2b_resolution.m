function result=audit_two_dimensional_closeout_s2b_resolution(out)
% Independent resolution classification and detailed point CSV, no new fields.
dest=fullfile(out,'formal_reference_v2'); assert(~isfile(fullfile(dest,'audit/resolution_audit.json')));
S=readtable(fullfile(dest,'sensitivity_fields.csv'),'TextType','string');
R=readtable(fullfile(dest,'reference.csv'),'TextType','string');
Q=readtable(fullfile(dest,'resolutions.csv'),'TextType','string');
L=readtable(fullfile(dest,'fields_line.csv'),'TextType','string'); F=readtable(fullfile(dest,'fields_core.csv'),'TextType','string');
P=readtable(fullfile(out,'data/sensitivity_fields.csv'),'TextType','string');
cases=["flat","weak_low_K","region_II_low_K","strong_height_low_K"];
maxdiff=0; all_decisions=true; recomputed=zeros(height(Q),1);
for j=1:height(Q)
    c=find(cases==Q.case_name(j));
    if Q.support_id(j)=="M99_z3_411", T=L; ix=(1:411)'; else, T=F; end
    a=getraw(T,Q.case_name(j),Q.candidate(j)); b=getraw(T,Q.case_name(j),Q.reference(j)); den=norm(b);
    if Q.support_id(j)=="P11-P33"
        coords=T(T.case_name==Q.case_name(j)&T.method==Q.reference(j),:);
        sample=S(S.case_index==c&S.variant=="FullFK_complete_GL800",:); ix=zeros(9,1);
        for k=1:9, ix(k)=sample.sample_index(sample.x_m==coords.x_m(k)&sample.z_m==coords.z_m(k)); end
    end
    components=zeros(1,4);
    if (Q.candidate(j)=="BIE"||Q.reference(j)=="BIE")&&c>1
        rr=R(R.case_name==Q.case_name(j),:); components(1)=norm(rr.resolution_estimate(ix))/den;
    end
    if Q.candidate(j)=="PE"||Q.reference(j)=="PE"
        pp=P(P.job_id=="pe"+string(c),:);
        u=varfield(pp,"primary"); ug=varfield(pp,"N1968_same_source"); half=varfield(pp,"half_step_line411");
        pe=getraw(L,Q.case_name(j),"PE"); components(2)=sqrt(numel(a))*(max(abs(ug-u))+max(abs(half-pe)))/den;
    end
    if Q.candidate(j)=="FullFK"||Q.reference(j)=="FullFK"
        s=S(S.case_index==c,:); hi=varfield(s,"FullFK_complete_GL800"); lo=varfield(s,"FullFK_complete_GL400");
        components(3)=norm(hi(ix)-lo(ix))/den;
    end
    flat=getraw(T,"flat","FullFK"); components(4)=256*eps*norm(flat)/den;
    recomputed(j)=sum(components); e=norm(a-b)/den;
    maxdiff=max([maxdiff,abs(recomputed(j)-Q.numerical_resolution_estimate(j)),abs(e-Q.relative_difference(j))]);
    all_decisions=all_decisions&&((e>5*recomputed(j))==logical(Q.resolved_above_5x_estimate(j)));
end
assert(maxdiff<1e-12 && all_decisions);
M=readtable(fullfile(dest,'point_metrics.csv'),'TextType','string');
candidate_abs=zeros(height(M),1); mag_sq=candidate_abs; phase_sq=candidate_abs; closure=candidate_abs;
phase_valid=false(height(M),1); near_null=phase_valid; near_null_empirical=phase_valid;
for j=1:height(M)
    a=getpoint(F,M.case_name(j),M.candidate(j),M.receiver_id(j)); b=getpoint(F,M.case_name(j),M.reference(j),M.receiver_id(j));
    candidate_abs(j)=abs(a); mag_sq(j)=(abs(a)-abs(b))^2; phase_sq(j)=4*abs(a)*abs(b)*sin(angle(a*conj(b))/2)^2;
    closure(j)=abs(a-b)^2-mag_sq(j)-phase_sq(j); phase_valid(j)=abs(a)>0&&abs(b)>0;
    support=getraw(F,M.case_name(j),M.reference(j)); near_null(j)=abs(b)<=1e-12*max(abs(support));
    if M.case_name(j)~="flat" && M.reference(j)=="BIE"
        rr=R(R.case_name==M.case_name(j),:); coord=F(F.case_name==M.case_name(j)&F.method==M.reference(j)&F.receiver_id==M.receiver_id(j),:);
        est=rr.resolution_estimate(rr.x_m==coord.x_m&rr.z_m==coord.z_m); near_null_empirical(j)=abs(b)<=5*est;
    end
end
detailed=addvars(M,candidate_abs,mag_sq,phase_sq,closure,phase_valid,near_null,near_null_empirical);
assert(max(abs(closure))<1e-12 && ~any(near_null) && ~any(near_null_empirical));
writetable(detailed,fullfile(dest,'point_metrics_detailed.csv')); save(fullfile(dest,'point_metrics_detailed.mat'),'detailed','-v7.3');
result=struct('pass',true,'resolution_rows',height(Q),'max_independent_resolution_error',maxdiff, ...
    'all_5x_classifications_match',all_decisions,'point_rows',height(detailed),'max_point_decomposition',max(abs(closure)), ...
    'near_null_scale_count',nnz(near_null),'near_null_empirical_count',nnz(near_null_empirical),'new_physical_fields',0);
closeout_s2b_json(fullfile(dest,'audit/resolution_audit.json'),result); disp(result);
end
function u=getraw(T,c,m)
r=T(T.case_name==c&T.method==m,:); u=complex(r.u_real,r.u_imag);
end
function u=getpoint(T,c,m,id)
r=T(T.case_name==c&T.method==m&T.receiver_id==id,:); u=complex(r.u_real,r.u_imag);
end
function u=varfield(T,name)
r=T(T.variant==name,:);
if ismember('sample_index',T.Properties.VariableNames), r=sortrows(r,'sample_index');
else, r=sortrows(r,'probe_index');
end
u=complex(r.u_real,r.u_imag);
end
