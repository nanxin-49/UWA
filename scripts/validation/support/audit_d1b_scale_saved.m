function receipt = audit_d1b_scale_saved(file,profile,src,workers,reference,monitor)
%AUDIT_D1B_SCALE_SAVED Independent GL112 and saved-system linear residual.
% Does not call a solver. The matrix is read in 32-row slabs from disk.
arguments
    file (1,1) string
    profile (1,1) struct
    src (1,1) struct
    workers (1,1) double
    reference (1,1) string
    monitor (1,1) struct
end
d1b_scale_phase(monitor,'saved_linear_audit'); t0=tic;
s=load(file,'b','rhs'); b=s.b; rhs=s.rhs; f=matfile(file);
assert(b.finite && all(isfinite(b.density),'all'));
assert(b.linear_residual<1e-10 && b.offgrid_boundary_residual<1e-4);
assert(isequal(b.receiver_x_m,(-1139:1139)'*(200/1139)) && isequal(b.receiver_z_m,6*ones(2279,1)));
residual=zeros(size(rhs),'like',rhs);
for first=1:32:b.node_count
    ii=first:min(first+31,b.node_count); block=f.A(ii,:);
    assert(all(isfinite(block),'all'));
    residual(ii,:)=block*b.density-rhs(ii,:);
end
linear=norm(residual)/max(norm(rhs),realmin); assert(linear<1e-10);
saved_linear_audit_s=toc(t0);
d1b_scale_phase(monitor,'independent_GL112'); t0=tic;
[ind,owner,ind_relative]=d1b_boundary_audit(b,profile,src,workers);
assert(all(isfinite(ind),'all') && isfinite(ind_relative) && ind_relative<1e-4);
independent_audit_s=toc(t0);
equivalence=struct('scope','1536 full serial solution; all base rows, fixed-density samples and all GL112 rows at larger sizes', ...
    'full_large_serial_solve_performed',false,'checks',[]);
d1b_scale_phase(monitor,'serial_equivalence');
if strlength(reference)>0
    r=load(reference,'b','rhs'); fi=matfile(reference);
    assert(b.node_count==1536 && r.b.node_count==1536);
    e2=0; scale2=0; peak=0; refpeak=0;
    for first=1:32:b.node_count
        ii=first:min(first+31,b.node_count); base=fi.A(ii,:); actual=f.A(ii,:);
        e2=e2+sum(abs(actual-base).^2,'all'); scale2=scale2+sum(abs(base).^2,'all');
        peak=max(peak,max(abs(actual-base),[],'all')); refpeak=max(refpeak,max(abs(base),[],'all'));
    end
    checks=struct('matrix_L2',sqrt(e2)/max(sqrt(scale2),realmin),'matrix_max',peak/max(refpeak,1));
    assert(checks.matrix_L2<1e-11 && checks.matrix_max<1e-11);
    checks.rhs=local_compare(rhs,r.rhs,1e-11);
    checks.density=local_compare(b.density,r.b.density,1e-10);
    checks.receiver=local_compare(b.receiver_field,r.b.receiver_field,1e-11);
    checks.midpoint=local_compare(b.boundary_total_residual,r.b.boundary_total_residual,1e-10);
    baseline=load(fullfile(fileparts(reference),'independent.mat'),'ind');
    checks.independent=local_compare(ind,baseline.ind,1e-10);
    equivalence.checks=checks;
end
assembly=b.serial_assembly_audit;
assert(assembly.relative_L2<1e-11 && assembly.max_scaled_abs<1e-11 && assembly.rows==b.node_count);
independent_worker_checks=struct('pid',{},'check',{});
if workers==2
    ids=sort(monitor.worker_ids(:));
    assert(isequal(sort([assembly.workers.pid].'),ids));
    for name={'assembly','receiver','boundary'}
        assert(isequal(sort(unique(b.worker_rows.(name{1}))),ids));
    end
    assert(isequal(sort(unique(owner)),ids));
    d1b_scale_phase(monitor,'serial_GL112_check');
    [serial_ind,~,serial_relative]=d1b_boundary_audit(b,profile,src,0);
    assert(serial_relative<1e-4);
    for w=1:workers
        ii=owner==ids(w);
        independent_worker_checks(w)=struct('pid',ids(w),'check',local_compare(ind(ii),serial_ind(ii),1e-10));
    end
end
save(fullfile(fileparts(file),'independent.mat'),'ind','owner','ind_relative','-v7');
receipt=struct('status','PASS_INDEPENDENT_SAVED_AUDIT','failed_checks',0, ...
    'nodes',b.node_count,'workers',workers,'linear_relative',linear, ...
    'midpoint_relative',b.offgrid_boundary_residual,'independent_relative',ind_relative, ...
    'saved_linear_audit_s',saved_linear_audit_s,'independent_audit_s',independent_audit_s, ...
    'serial_assembly_audit',assembly,'equivalence',equivalence,'timings_s',b.cost_meta.timings, ...
    'fixed_density_serial_checks',b.fixed_density_serial_checks,'GL112_serial_worker_checks',independent_worker_checks, ...
    'data_sha256',d1b_scale_hash(file));
d1b_scale_json(fullfile(fileparts(file),'audit.json'),receipt);
end

function v=local_compare(a,b,tol)
assert(isequal(size(a),size(b)) && all(isfinite(a),'all') && all(isfinite(b),'all'));
v=struct('relative_L2',norm(a(:)-b(:))/max(norm(b(:)),realmin), ...
    'max_scaled_abs',max(abs(a-b),[],'all')/max(max(abs(b),[],'all'),1));
assert(v.relative_L2<tol && v.max_scaled_abs<tol);
end
