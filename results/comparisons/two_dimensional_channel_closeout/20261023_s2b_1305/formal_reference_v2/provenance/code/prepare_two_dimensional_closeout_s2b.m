function plan=prepare_two_dimensional_closeout_s2b(out,available_physical_bytes)
%PREPARE_TWO_DIMENSIONAL_CLOSEOUT_S2B Read-only P1 and frozen job contract.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(isfile(fullfile(out,'data/p0_audit.json')) && ~isfile(fullfile(out,'data/p1.mat')));
assert(jsondecode(fileread(fullfile(out,'data/p0_audit.json'))).pass);
old=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final');
t=load(fullfile(out,'data/p0_data.mat'),'fields'); fields=t.fields;
manifest=readtable(fullfile(old,'data/source_manifest.csv'),'TextType','string');
for j=1:height(manifest)
    assert(strcmpi(closeout_s2b_hash(manifest.path(j)),manifest.sha256(j)),'S2a input/code drift');
end
base=fullfile(root,'results/visualization/pe_bellhop_bie_xz/20261019_115214_348/source_data');
files={fullfile(base,'weak_reference.mat'), ...
    fullfile(root,'results/archive/cl-verification-worktree-unmanaged-20260930/validation/pe_bellhop_helmholtz_bie_reference/R4_region_II/R4_region_II_validation.mat'), ...
    fullfile(base,'strong_reference.mat')};
snapshots=cell(3,1); rows={}; amplitudes=[.01 .05 .2]; counts=[3 3 4];
for c=1:3
    assert(isfile(files{c}));
    names={'frequency_hz','c0_mps','nx','xw_m','sigma_src_m','z_tx_m','surface_support_m','surface_taper_inner_m','amplitude_m'};
    expected=[4000 1500 984 192.1875 .3 100 50 42 amplitudes(c)];
    for j=1:numel(names), assert(h5read(files{c},['/validation/config/' names{j}])==expected(j)); end
    suites={'spatial_runs','window_runs'}; sizes=[counts(c) 3];
    if c==1, suites{3}='quadrature_runs'; sizes(3)=2; end
    ss=struct;
    for s=1:numel(suites)
        levels=cell(sizes(s),1);
        for j=1:sizes(s)
            b=closeout_read_bie_snapshot(files{c},suites{s},j);
            [eta,ep]=closeout_s2b_surface(b.node_x_m,amplitudes(c));
            normal=[-ep./hypot(1,ep),1./hypot(1,ep)];
            err=max(abs(eta-b.node_z_m)); ne=max(abs(normal-b.node_normal),[],'all');
            assert(err<1e-12 && ne<1e-12 && b.k_radpm==2*pi*4000/1500);
            assert(all(isfinite(b.density)) && b.halfplane_h_m<min(b.node_z_m));
            levels{j}=b;
            rows{end+1,1}=table(c,string(suites{s}),j,b.node_count,b.half_width_m, ...
                b.points_per_wavelength_actual,b.self_quadrature_order,b.linear_residual,b.offgrid_boundary_residual, ...
                err,ne,string(files{c}),string(b.hdf5_group), ...
                'VariableNames',{'case_index','suite','level','node_count','half_width_m','ppw','self_quad', ...
                'linear_residual','boundary_residual','geometry_error','normal_error','file','hdf5_group'}); %#ok<AGROW>
        end
        ss.(suites{s})=levels;
    end
    fine=ss.spatial_runs{end}; [ok,ix]=ismember(fields.x_m,fine.receiver_x_m);
    assert(all(ok)); saved=fine.receiver_field(ix);
    assert(norm(saved-fields.BIE{c+1})/norm(saved)<1e-11,'BIE saved line mismatch');
    snapshots{c}=ss;
end
% Charge every persisted S2a run, including diagnostic runs; do not trust prose alone.
prior=0; costrows={};
for nm={'20261023_s2a_a','20261023_s2a_b','20261023_s2a_final'}
    f=fullfile(root,'results/comparisons/two_dimensional_channel_closeout',nm{1},'data/cost.json');
    a=jsondecode(fileread(f)); sec=a.plane_seconds+a.PE_flat_seconds+sum(a.PE_rough_seconds)+sum(a.FullFK_seconds);
    prior=prior+sec; costrows{end+1,1}=table(string(nm{1}),sec,'VariableNames',{'run_id','recorded_core_seconds'}); %#ok<AGROW>
end
% Untimed historical audit integrations/pre-output attempt: reserve, not measured cost.
unmetered_prior_reserve=120;
estimated_density_peak=64*1024^2; estimated_dense_resolve_peak=4*16*4480^2+128*1024^2;
assert(estimated_density_peak<=.7*available_physical_bytes);
plan=struct('status','P1_INPUT_RESOURCE_GATE_PASS_DATA_EVALUATION_ONLY', ...
    'core_limit_seconds',7200,'single_limit_seconds',900,'reserve_fraction',.30, ...
    'prior_recorded_seconds',prior,'prior_unmetered_reserve_seconds',unmetered_prior_reserve, ...
    'remaining_seconds',7200-prior-unmetered_prior_reserve, ...
    'available_physical_bytes',available_physical_bytes,'density_peak_estimate_bytes',estimated_density_peak, ...
    'dense_resolve_peak_estimate_bytes',estimated_dense_resolve_peak, ...
    'dense_resolve_ram_gate',estimated_dense_resolve_peak<=.7*available_physical_bytes, ...
    'new_boundary_solves_authorized',false,'density_evaluation_preestimate_seconds',5, ...
    'dense_resolve_seconds','unmeasured; not authorized without a separate representative solve gate', ...
    'reference_status','REFERENCE_INCONCLUSIVE_until_matched_point_sensitivity', ...
    'version_limit','Saved formulation/schema and numeric config verified; original execution-time solver source hash absent.', ...
    'scope','S2b only; no S3-S6/Bellhop/communication', ...
    'reference_rule','max absolute field changes across spatial/window/available quadrature levels; empirical sensitivity, not rigorous bound; separation >=5*sensitivity only diagnostic', ...
    'numerical_jobs','BIE all saved levels at fixed 9 points and line extrema; PE half-step and doubled transverse sampling with same frozen spectrum; FullFK 2049/4097/8193 at L60 and L80 same spacing; no new density solve');
coverage=readtable(fullfile(old,'data/s2b_coverage_gap.csv'),'TextType','string');
assert(nnz(coverage.available)==48); inventory=vertcat(rows{:}); prior_cost=vertcat(costrows{:});
save(fullfile(out,'data/p1.mat'),'plan','snapshots','coverage','inventory','prior_cost','-v7.3');
writetable(inventory,fullfile(out,'data/bie_density_inventory.csv'));
writetable(prior_cost,fullfile(out,'data/prior_costs.csv'));
writetable(coverage,fullfile(out,'data/coverage_before.csv'));
closeout_s2b_json(fullfile(out,'data/p1_gate.json'),plan);
extra=table(string(files(:)),strings(3,1),'VariableNames',{'path','sha256'});
for j=1:3, extra.sha256(j)=string(closeout_s2b_hash(files{j})); end
writetable([manifest;extra],fullfile(out,'data/input_hashes.csv'));
disp(plan);
end
