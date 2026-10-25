function audit = audit_two_dimensional_closeout_s2c_b_phi_final(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_PHI_FINAL Independent PHI saved-CSV audit.
% No propagation, boundary solve/evaluation, or Kirchhoff integration.
arguments
    run_dir (1,:) char
end
assert(isfolder(run_dir) && contains(run_dir,'20261023_s2c_b_phi_plus90_01'));
analysis=fullfile(run_dir,'analysis'); dest=fullfile(run_dir,'audit');
target=fullfile(dest,'final_independent_audit_phi_v2.json');
assert(~isfile(target),'Independent audit already exists; do not overwrite');
R=readtable(fullfile(analysis,'raw_fields.csv'));
P=readtable(fullfile(analysis,'model_point_diagnostics.csv'));
S=readtable(fullfile(analysis,'support_metrics.csv'));
Q=readtable(fullfile(analysis,'paired_changes.csv'));
H=readtable(fullfile(analysis,'BIE_physical_pair.csv'));
V=readtable(fullfile(run_dir,'P1','BIE_reference_versions.csv'));
B=readtable(fullfile(run_dir,'P1','BIE_resolution_point.csv'));
assert(height(R)==2502 && height(P)==1668 && height(S)==20 && height(Q)==10);
assert(height(H)==417 && height(V)==2085 && height(B)==417);
cases=["BASE_M","PHI_PLUS90"]; methods=["PE","BIE","FullFK"];
U=complex(zeros(417,2,3)); X=complex(zeros(417,2,3));
for c=1:2
    for m=1:3
        a=R(string(R.case_name)==cases(c) & string(R.method)==methods(m),:);
        assert(height(a)==417 && isequal(a.sample_index,(1:417)'));
        if c==1 && m==1
            x=a.x_m(:); z=a.z_m(:);
        else
            assert(isequal(a.x_m(:),x) && isequal(a.z_m(:),z));
        end
        assert(all(a.valid) && all(isfinite(a.u_native_real)) && all(isfinite(a.u_native_imag)));
        U(:,c,m)=complex(a.u_native_real,a.u_native_imag);
        X(:,c,m)=complex(a.u_flat_real,a.u_flat_imag);
        assert(max(abs(complex(a.delta_u_real,a.delta_u_imag)-(U(:,c,m)-X(:,c,m))))<2e-14);
        assert(max(abs(complex(a.G_real,a.G_imag)-U(:,c,m)./X(:,c,m)))<2e-13);
    end
end
assert(all(z(1:411)==3) && size(unique([x,z],'rows'),1)==417);
assert(isequal(H.x_m(:),x) && isequal(H.z_m(:),z));
assert(isequal(B.x_m(:),x) && isequal(B.z_m(:),z));
assert(all(abs(x(1:411))<=41));
assert(max(abs(U(:,1,2)-complex(H.BASE_BIE_real,H.BASE_BIE_imag)))<2e-14);
assert(max(abs(U(:,2,2)-complex(H.PHI_BIE_real,H.PHI_BIE_imag)))<2e-14);
assert(max(abs(complex(H.delta_phi_real,H.delta_phi_imag) ...
    -(U(:,2,2)-U(:,1,2))))<2e-14);
assert(max(abs(H.delta_phi_abs-abs(U(:,2,2)-U(:,1,2))))<2e-14);
v=complex(zeros(417,5)); ids=["main96","self64","space14","window60","window70"];
for j=1:5
    w=V(string(V.level)==ids(j),:);
    assert(height(w)==417 && isequal(w.sample_index,(1:417)'));
    assert(isequal(w.x_m(:),x) && isequal(w.z_m(:),z));
    v(:,j)=complex(w.u_real,w.u_imag);
end
bres=abs(v(:,1)-v(:,2))+abs(v(:,1)-v(:,3))+abs(v(:,4)-v(:,5));
assert(max(abs(v(:,1)-U(:,2,2)))<2e-14);
assert(max(abs(bres-B.empirical_resolution_abs))<2e-14);
assert(max(abs(H.PHI_BIE_resolution_abs-bres))<2e-14);
ix9=find(abs(x)<=.2 & abs(z-3)<=.126);
assert(numel(ix9)==9 && numel(unique(x(ix9)))==3 && numel(unique(z(ix9)))==3);
center=find(abs(x(1:411))<=1); edge=find(abs(x(1:411))>=30);
cross=ix9(x(ix9)==0 | z(ix9)==3);
ii={(1:411)',center,ix9,cross,edge};
names=["line_411","center_11","nine_9","cross_5","edge_104"];
assert(isequal(cellfun(@numel,ii),[411 11 9 5 104]));
point_error=0; identity_error=0; support_error=0; pair_error=0;
for c=1:2
    ref=U(:,c,2);
    for m=[1 3]
        u=U(:,c,m); e=u-ref; ph=angle(u.*conj(ref));
        amp2=(abs(u)-abs(ref)).^2;
        phase2=4*abs(u).*abs(ref).*sin(ph/2).^2;
        p=P(string(P.case_name)==cases(c) & string(P.method)==methods(m),:);
        assert(height(p)==417 && isequal(p.sample_index,(1:417)'));
        point_error=max([point_error;abs(complex(p.residual_real,p.residual_imag)-e); ...
            abs(p.residual_abs-abs(e));abs(p.amplitude_difference_dB-20*log10(abs(u)./abs(ref))); ...
            abs(p.wrapped_phase_difference_rad-ph);abs(p.amplitude_square-amp2); ...
            abs(p.phase_square-phase2);abs(p.complex_square-abs(e).^2)]);
        identity_error=max(identity_error,max(abs(abs(e).^2-amp2-phase2)));
        hb=[H.BASE_BIE_resolution_abs,H.PHI_BIE_resolution_abs];
        assert(max(abs(p.BIE_resolution_abs-hb(:,c)))<2e-14);
        for k=1:5
            t=S(string(S.case_name)==cases(c)&string(S.method)==methods(m) ...
                &string(S.support_id)==names(k),:);
            assert(height(t)==1 && t.n==numel(ii{k}));
            idx=ii{k}; den=mean(abs(ref(idx)).^2);
            L=sqrt(mean(abs(e(idx)).^2)/den);
            A=sqrt(mean(amp2(idx))/den); F=sqrt(mean(phase2(idx))/den);
            vals=[den,L,A,F,mean(abs(e(idx))),median(abs(e(idx))),max(abs(e(idx)))];
            old=[t.BIE_denominator,t.complex_L2,t.amplitude_L2,t.phase_L2, ...
                t.mean_abs_residual,t.median_abs_residual,t.max_abs_residual];
            support_error=max(support_error,max(abs(vals-old)));
            assert(t.near_zero_count==sum(p.reference_near_zero(idx)));
            hb=[H.BASE_BIE_resolution_abs,H.PHI_BIE_resolution_abs];
            assert(abs(t.BIE_resolution-norm(hb(idx,c))/norm(ref(idx)))<2e-12);
            assert(t.resolved_above_5x==(L>5*t.total_empirical_resolution));
        end
    end
end
for k=1:5
    idx=ii{k}; b0=U(idx,1,2); b1=U(idx,2,2);
    phys=norm(b1-b0)/norm(b0);
    for m=[1 3]
        t=Q(string(Q.method)==methods(m)&string(Q.support_id)==names(k),:);
        assert(height(t)==1 && t.n==numel(idx));
        e0=U(idx,1,m)-b0; e1=U(idx,2,m)-b1;
        vals=[mean(abs(b0).^2),phys,norm(e0)/norm(b0), ...
            norm(e1)/norm(b1),norm(e1)/norm(b1)-norm(e0)/norm(b0), ...
            norm(e1-e0)/norm(b0)];
        old=[t.BASE_BIE_denominator,t.BIE_physical_response_L2, ...
            t.BASE_model_error_L2,t.PHI_model_error_L2, ...
            t.signed_model_error_change,t.paired_residual_change_L2];
        pair_error=max(pair_error,max(abs(vals-old)));
        assert(t.BIE_physical_resolved_above_5x==(phys>5*t.BIE_paired_resolution));
        assert(t.paired_residual_change_resolved_above_5x== ...
            (vals(6)>5*t.paired_residual_resolution));
    end
end
assert(point_error<2e-12 && identity_error<2e-12);
assert(support_error<2e-12 && pair_error<2e-12);
% Sensitivity provenance: PHI_PLUS90 five BIE levels, PE half step/probes,
% and FullFK GL8/extent, directly from saved numerical jobs.
J=load(fullfile(run_dir,'jobs','pe_phi.mat'),'payload','receipt');
F=load(fullfile(run_dir,'jobs','fullfk_phi.mat'),'payload','receipt');
assert(strcmp(F.payload.version,'FullFK_complete_analytic_flat_exterior_GL8_v1'));
assert(J.receipt.flat_control_relative<1e-10);
assert(strcmp(J.payload.surface_version, ...
    'C2 PHI_PLUS90 A0.05 K0.1 phi+pi/2 physical depth'));
assert(contains(J.payload.double_grid_source,'zero-padded SAME frozen N984'));
assert(F.payload.GL8_800.local_radius_m==Inf);
assert(F.payload.component_max_abs<1e-12);
assert(max(abs(J.payload.u_native(:)-U(:,2,1)))<2e-14);
assert(max(abs(F.payload.u_native(:)-U(:,2,3)))<2e-14);
pe_step=max(abs(J.payload.half_step_native(:)-J.payload.u_native(:)));
pe_probe=max(abs(J.payload.double_grid_probe_native(:)-J.payload.primary_probe_native(:)));
fk_quad=norm(F.payload.GL8_800.receiver_field(:)-F.payload.GL8_400.receiver_field(:))/norm(U(:,2,3));
fk_extent=norm(F.payload.extent60.receiver_field(:)-F.payload.extent90.receiver_field(:))/norm(U(:,2,3));
assert(pe_step>0 && pe_probe>0 && fk_quad<1e-10 && fk_extent<1e-10);
% Independently rebuild both-case empirical resolution components and all
% five support classifications from the saved numerical runs.
formal=load(fullfile(root_of_run(run_dir),'results','comparisons', ...
    'two_dimensional_channel_closeout','20261023_s2b_1305', ...
    'formal_reference_v2','reference_data.mat'),'fk_runs','numerical','flat_exact');
old=load(fullfile(root_of_run(run_dir),'results','comparisons', ...
    'two_dimensional_channel_closeout','20261023_s2b_1305', ...
    'checkpoints','pe3.mat'),'payload');
base_ref=H.BASE_BIE_resolution_abs; phi_ref=bres;
base_pe=max(abs(old.payload.step_half_line(:)-U(1:411,1,1))) ...
    +max(abs(old.payload.grid_double_u(:)-old.payload.u(:)));
phi_pe=pe_step+pe_probe;
base_fk=abs(formal.fk_runs{3,2}.receiver_field(:) ...
    -formal.fk_runs{3,1}.receiver_field(:));
tail=formal.numerical.FullFK_independent_tail_max_abs( ...
    formal.numerical.case_name=="region_II_low_K");
assert(isscalar(tail)); base_fk=base_fk+tail;
phi_fk=abs(F.payload.GL8_800.receiver_field(:)-F.payload.GL8_400.receiver_field(:)) ...
    +abs(F.payload.extent60.receiver_field(:)-F.payload.extent90.receiver_field(:)) ...
    +abs(F.payload.GL8_800.receiver_field(:)-F.payload.extent60.receiver_field(:));
br=[base_ref,phi_ref]; fr=[base_fk,phi_fk]; pr=[base_pe,phi_pe];
resolution_error=0; classification_count=0;
for k=1:5
    idx=ii{k}; n=numel(idx); den0=norm(U(idx,1,2));
    phys_floor=(norm(br(idx,1))+norm(br(idx,2)))/den0;
    for m=[1 3]
        total=zeros(1,2);
        for c=1:2
            den=norm(U(idx,c,2));
            bb=norm(br(idx,c))/den;
            pp=(m==1)*sqrt(n)*pr(c)/den;
            ff=(m==3)*norm(fr(idx,c))/den;
            rr=256*eps*norm(formal.flat_exact(idx))/den;
            total(c)=bb+pp+ff+rr;
            t=S(string(S.case_name)==cases(c)&string(S.method)==methods(m) ...
                &string(S.support_id)==names(k),:);
            stored=[t.BIE_resolution,t.PE_probe_resolution, ...
                t.FullFK_resolution,t.roundoff_allowance,t.total_empirical_resolution];
            resolution_error=max(resolution_error,max(abs(stored-[bb,pp,ff,rr,total(c)])));
            assert(t.resolved_above_5x==(t.complex_L2>5*total(c)));
            classification_count=classification_count+1;
        end
        t=Q(string(Q.method)==methods(m)&string(Q.support_id)==names(k),:);
        pair_floor=(norm(br(idx,1))+norm(br(idx,2)) ...
            +(m==1)*sqrt(n)*sum(pr) ...
            +(m==3)*(norm(fr(idx,1))+norm(fr(idx,2))) ...
            +2*256*eps*norm(formal.flat_exact(idx)))/den0;
        resolution_error=max(resolution_error,max(abs([t.BIE_paired_resolution-phys_floor, ...
            t.paired_residual_resolution-pair_floor])));
        assert(t.BIE_physical_resolved_above_5x==(t.BIE_physical_response_L2>5*phys_floor));
        assert(t.scalar_error_change_resolved_above_5x== ...
            (abs(t.signed_model_error_change)>5*sum(total)));
        assert(t.paired_residual_change_resolved_above_5x== ...
            (t.paired_residual_change_L2>5*pair_floor));
        classification_count=classification_count+3;
    end
end
assert(resolution_error<2e-12 && classification_count==50);
jobs=["main96","self64","space14","window60","window70","pe_phi","fullfk_phi"];
cost=zeros(numel(jobs),1); ram=zeros(numel(jobs),1);
src_sha=''; cfg_sha=''; solver_sha='';
gate=closeout_s2c_b_phi_gate(run_dir);
for j=1:numel(jobs)
    if j<=5, name="bie_"+jobs(j); else, name=jobs(j); end
    receipt=jsondecode(fileread(fullfile(run_dir,'jobs',char(name+".json"))));
    assert(strcmp(receipt.status,'COMPLETE') && receipt.wall_seconds<=900);
    assert(receipt.estimated_increment_bytes<=.7*receipt.available_physical_bytes);
    assert(receipt.authorized_subbudget_s==gate.subbudget_s);
    assert(strcmp(receipt.budget_amendment_sha256,gate.budget_amendment_sha256));
    if j<=5
        if j==1
            src_sha=receipt.source_sha256; cfg_sha=receipt.configuration_sha256;
            solver_sha=receipt.solver_sha256;
        else
            assert(strcmp(src_sha,receipt.source_sha256));
            assert(strcmp(cfg_sha,receipt.configuration_sha256));
            assert(strcmp(solver_sha,receipt.solver_sha256));
        end
        assert(receipt.receiver_values==417 && receipt.linear_residual<1e-8);
    end
    cost(j)=receipt.wall_seconds; ram(j)=receipt.available_physical_bytes;
end
root=root_of_run(run_dir);
assert(strcmp(src_sha,'8d60d3d9d6aa3aa16110ac33d81e0d193a61526aef2c0c2e9edd531cae531f55'));
assert(strcmp(cfg_sha,gate.configuration_sha256));
assert(strcmp(solver_sha,closeout_s2b_hash(fullfile(root,'scripts','validation', ...
    'support','solve_helmholtz_bie_halfplane_vertical.m'))));
ref_audit=jsondecode(fileread(fullfile(dest,'P1_five_BIE_reference_audit.json')));
assert(strcmp(ref_audit.status,'FIVE_BIE_LEVELS_COMPLETE_EMPIRICAL_REFERENCE'));
assert(strcmp(ref_audit.case_id,'PHI_PLUS90') && ref_audit.BIE_point_count==417);
assert(sum(cost)<=gate.subbudget_s && ...
    gate.shared_cap_s-gate.charged_prior_seconds-sum(cost)>=gate.S3A_reserve_s);
% Redraw from saved MAT through a fresh reporting invocation, then compare
% exact PNG pixels with the delivered figures. FIG metadata is reopened.
addpath(fullfile(root,'scripts','reporting'));
redraw=fullfile(dest,'independent_redraw_phi_v2');
render_two_dimensional_closeout_s2c_b_phi(fullfile(analysis,'analysis_data.mat'),redraw);
geometry_redraw=fullfile(dest,'independent_geometry_redraw_phi_v2');
render_two_dimensional_closeout_s2c_b_phi_geometry(root,geometry_redraw, ...
    fullfile(run_dir,'figures_final','geometry.mat'));
stems=["M05a_geometry_only","M05b_BIE_physical_response", ...
    "M05c_model_error_pair","M05d_line_center_nine","M05e_mirror_symmetry"];
axes_min=[3 9 8 8 6];
pixel_ok=false(5,1); fig_ok=false(5,1); dims=zeros(5,2);
for j=1:5
    name=char(stems(j)); original=fullfile(run_dir,'figures_final_complete_v2',[name '.png']);
    image_dir=redraw; if j==1, image_dir=geometry_redraw; end
    replica=fullfile(image_dir,[name '.png']);
    a=imread(original); b=imread(replica);
    dims(j,:)=size(a,[1 2]); pixel_ok(j)=isequal(a,b);
    f=openfig(fullfile(image_dir,[name '.fig']),'invisible');
    axes_count=numel(findall(f,'Type','axes'));
    fig_ok(j)=axes_count>=axes_min(j) && numel(findall(f,'Type','text'))>0;
    close(f);
end
assert(all(pixel_ok) && all(fig_ok));
M=readtable(fullfile(run_dir,'figures_final_complete_v2','M05e_mirror_pairs.csv'));
assert(height(M)==2502);
[paired,mirror]=ismember([-x,z],[x,z],'rows'); assert(all(paired));
mirror_error=0;
for c=1:2
    for m=1:3
        q=M(string(M.case_name)==cases(c)&string(M.method)==methods(m),:);
        assert(height(q)==417 && isequal(q.sample_index(:),(1:417)'));
        assert(isequal(q.x_m(:),x) && isequal(q.z_m(:),z));
        assert(isequal(q.mirror_index(:),mirror));
        u=U(:,c,m); e=u-U(:,c,2);
        stored=complex(q.native_mirror_difference_real,q.native_mirror_difference_imag);
        residual=complex(q.residual_mirror_difference_real,q.residual_mirror_difference_imag);
        mirror_error=max([mirror_error; ...
            abs(complex(q.u_real,q.u_imag)-u); ...
            abs(complex(q.mirror_u_real,q.mirror_u_imag)-u(mirror)); ...
            abs(stored-(u-u(mirror))); ...
            abs(residual-(e-e(mirror)))]);
    end
end
assert(mirror_error<2e-14);
audit=struct('status','S2C_B_PHI_PLUS90_INDEPENDENT_AUDIT_PASS', ...
    'raw_rows',height(R),'point_rows',height(P),'support_rows',height(S), ...
    'paired_rows',height(Q),'BIE_version_rows',height(V), ...
    'max_point_error',point_error,'max_identity_error',identity_error, ...
    'max_support_error',support_error,'max_paired_error',pair_error, ...
    'max_resolution_error',resolution_error,'classification_checks',classification_count, ...
    'PHI_PE_step_max_abs',pe_step,'PHI_PE_frozen_probe_max_abs',pe_probe, ...
    'PHI_FullFK_quadrature_relative',fk_quad,'PHI_FullFK_extent_relative',fk_extent, ...
    'job_seconds',cost,'new_acoustic_seconds',sum(cost), ...
    'minimum_job_available_physical_bytes',min(ram), ...
    'charged_core_seconds',gate.charged_prior_seconds+sum(cost), ...
    'S3A_planning_balance_seconds',gate.shared_cap_s-gate.charged_prior_seconds-sum(cost), ...
    'redraw_pixel_equal',pixel_ok,'fig_reopened',fig_ok,'image_dimensions_px',dims, ...
    'mirror_pair_rows',height(M),'mirror_csv_max_error',mirror_error, ...
    'source_sha256',src_sha,'configuration_sha256',cfg_sha, ...
    'BIE_solver_sha256',solver_sha, ...
    'qualification','empirical five-level BIE; PE grid sensitivity only 13 frozen probes');
fid=fopen(target,'w'); assert(fid>0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(audit,PrettyPrint=true));
disp(audit);
end

function root=root_of_run(run_dir)
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
end
