function audit = audit_two_dimensional_closeout_s2c_b_p0(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_P0 Independent zero-physics audit.
arguments
    run_dir (1,:) char
end
assert(isfolder(run_dir)); dest=fullfile(run_dir,'audit'); assert(~isfolder(dest));
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
formal=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2');
prior=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642');
t=load(fullfile(run_dir,'P0','p0_data.mat'));
g=readtable(fullfile(run_dir,'P0','geometry.csv'),TextType='string');
assert(height(g)==2*numel(t.gx));
a=g(g.case_id=="BASE_M",:); b=g(g.case_id=="K_HIGH",:);
geo_delta=max(abs([a.eta_m-t.eta0;a.slope-t.d10;a.curvature_per_m-t.d20; ...
    b.eta_m-t.eta1;b.slope-t.d11;b.curvature_per_m-t.d21]));
assert(geo_delta<1e-12);
assert(all(b.eta_m(abs(b.x_m)>=50)==0) && all(b.slope(abs(b.x_m)>=50)==0) ...
    && all(b.curvature_per_m(abs(b.x_m)>=50)==0));
raw=readtable(fullfile(run_dir,'P0','BASE_M_fields.csv'));
assert(height(raw)==417);
d=load(fullfile(formal,'reference_data.mat'),'U','x','z');
base=squeeze(d.U(:,3,:));
recovered=[complex(raw.PE_real,raw.PE_imag),complex(raw.BIE_real,raw.BIE_imag), ...
    complex(raw.FullFK_real,raw.FullFK_imag)];
field_delta=max(abs(recovered-base),[],'all');
assert(field_delta<1e-14 && max(abs(t.Ubase-base),[],'all')==0);
assert(max(abs(raw.x_m-d.x))==0 && max(abs(raw.z_m-d.z))==0);
cfg=fullfile(prior,'protocol','K_HIGH_configuration.json');
h=readtable(fullfile(prior,'protocol','configuration_hashes.csv'),TextType='string');
assert(strcmp(closeout_s2b_hash(cfg),h.sha256(h.case_id=="K_HIGH")));
authorization=jsondecode(fileread(fullfile(run_dir,'protocol','authorization.json')));
assert(strcmp(authorization.authorized_case,'K_HIGH') && ~authorization.other_candidates_authorized);
assert(strcmp(authorization.frozen_configuration_sha256,closeout_s2b_hash(cfg)));
gate=jsondecode(fileread(fullfile(run_dir,'protocol','P1_initial_gate.json')));
recalc=1.2*sum(gate.historical_stage_basis_s)+300;
assert(abs(recalc-gate.total_conservative_s)<1e-10);
assert(gate.new_acoustic_seconds==0);
assert(gate.RAM_pass == (gate.max_increment_bytes<=.7*gate.physical_available_bytes));
assert(gate.cost_pass && gate.total_conservative_s<=1800);
mkdir(dest);
render_two_dimensional_closeout_s2c_b_geometry(fullfile(run_dir,'P0','p0_data.mat'),fullfile(dest,'redraw'));
expected=imread(fullfile(run_dir,'figures_v2','M04a_geometry_only.png'));
actual=imread(fullfile(dest,'redraw','M04a_geometry_only.png'));
pixel_identical=isequal(expected,actual); assert(pixel_identical);
fig=openfig(fullfile(run_dir,'figures_v2','M04a_geometry_only.fig'),'invisible');
fig_reopen=~isempty(findall(fig,'Type','Axes')); close(fig); assert(fig_reopen);
audit=struct('pass',true,'P0_only',true,'new_acoustic_fields',0, ...
    'BASE_M_csv_mat_max_abs',field_delta,'geometry_csv_mat_max_abs',geo_delta, ...
    'K_HIGH_strictly_flat_outside_50',true,'configuration_hash_pass',true, ...
    'cost_recomputation_pass',true,'initial_RAM_pass',gate.RAM_pass, ...
    'figure_png_pixel_identical',pixel_identical,'figure_fig_reopened',fig_reopen, ...
    'source_S2b_status','S2B_COMPLETE_REFERENCE_QUALIFIED', ...
    'source_S2c_A_status','S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE');
closeout_s2b_json(fullfile(dest,'P0_independent_audit.json'),audit); disp(audit);
end
