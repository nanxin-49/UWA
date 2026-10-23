function receipt=run_two_dimensional_closeout_s2b_bie_reference(out,id,available_physical_bytes)
% Versioned BIE reference/cost work; invokes unchanged existing solver.
assert(any(strcmp(id,{'pilot2','pilot4','medium96','medium64','medium128','strong96','strong64','strong128','weak64','weak128'})));
dest=fullfile(out,'reference_followup_01'); assert(isfolder(dest));
file=fullfile(dest,['bie_' id '.mat']); assert(~isfile(file),'No duplicate solve');
p=load(fullfile(out,'data/p0_data.mat'),'fields'); src=p.fields.source;
main=load(fullfile(out,'data/s2b_data.mat'),'grid'); g=main.grid;
c=3; ppw=12; qa=128;
if strcmp(id,'pilot2'), ppw=2; qa=96;
elseif strcmp(id,'pilot4'), ppw=4; qa=96;
elseif strcmp(id,'medium96'), qa=96;
elseif strcmp(id,'medium64'), qa=64;
elseif strcmp(id,'strong96'), c=4; qa=96; ppw=14;
elseif strcmp(id,'strong64'), c=4; qa=64; ppw=14;
elseif strcmp(id,'strong128'), c=4;
elseif strcmp(id,'weak64'), c=2; qa=64;
elseif strcmp(id,'weak128'), c=2;
end
A=[0 .01 .05 .2]; amp=A(c); n=8*ceil(120/(.375/ppw*8));
estimated_peak=4*16*n^2+128*1024^2;
assert(estimated_peak<=.7*available_physical_bytes,'New BIE physical RAM gate failed');
if startsWith(id,'pilot')
    estimated_seconds=180;
else
    a=load(fullfile(dest,'bie_pilot2.mat'),'receipt'); b=load(fullfile(dest,'bie_pilot4.mat'),'receipt');
    % Observed whole-solve scaling plus independently measured dense solve.
    estimate_quad=max(a.receipt.wall_seconds*(n/a.receipt.node_count)^2, ...
        b.receipt.wall_seconds*(n/b.receipt.node_count)^2)*max(1,(qa/96)^3);
    lu_bound=b.receipt.linear_solve_micro_seconds*(n/1280)^3;
    estimated_seconds=2*(estimate_quad+lu_bound);
    if c==4
        % Refine the cost model only after the actual full 3840-node job.
        fullpilot=load(fullfile(dest,'bie_medium96.mat'),'receipt');
        estimated_seconds=2*(fullpilot.receipt.wall_seconds*(n/fullpilot.receipt.node_count)^2*max(1,(qa/96)^3)+lu_bound);
    end
    assert(estimated_seconds<900,'New BIE time gate failed');
end
cost=jsondecode(fileread(fullfile(out,'data/cost.json'))); charged=cost.charged_total_seconds;
prior=dir(fullfile(dest,'bie_*.mat'));
for j=1:numel(prior), q=load(fullfile(dest,prior(j).name),'receipt'); charged=charged+q.receipt.wall_seconds; end
assert(charged+estimated_seconds<7200,'Cumulative budget gate failed');
if ~startsWith(id,'pilot')
    probes=load(fullfile(out,'checkpoints',sprintf('bie%d.mat',c-1)),'payload');
    g=struct('x_m',probes.payload.x_m,'z_m',probes.payload.z_m);
end
cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',60, ...
    'points_per_wavelength',ppw,'panel_order',8,'self_quadrature_order',qa, ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
    'eta_fn',@(x) closeout_s2b_surface(x,amp), ...
    'eta_prime_fn',@(x) surface_derivative(x,amp), ...
    'incident_fn',@(x,z) fk_single_receiver_incident(src,x,z), ...
    'receiver_x_m',g.x_m,'receiver_z_m',g.z_m);
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
solver=fullfile(root,'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m');
solver_sha256=closeout_s2b_hash(solver);
% The original solver assumes vector receiver_z_m and handles it elementwise.
started=tic; boundary=solve_helmholtz_bie_halfplane_vertical(cfg); seconds=toc(started);
assert(seconds<900,'900 s cap exceeded: do not start more jobs');
cfg_numeric=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn'});
micro_seconds=NaN;
if strcmp(id,'pilot4')
    % Deterministic non-acoustic LU timing, no RNG state change.
    mm=1280; ii=(1:mm)'; bb=complex(1./(1+abs(ii-ii.')),1e-3*sin(ii+ii.')); bb=bb+mm*eye(mm);
    tt=tic; sol=bb\ones(mm,1); micro_seconds=toc(tt); assert(all(isfinite(sol))); clear bb sol
end
receipt=struct('job_id',id,'case_index',c,'amplitude_m',amp,'node_count',n,'ppw',ppw,'self_quadrature_order',qa, ...
    'wall_seconds',seconds,'linear_solve_micro_seconds',micro_seconds,'estimated_seconds',estimated_seconds, ...
    'available_physical_bytes',available_physical_bytes,'estimated_peak_increment_bytes',estimated_peak, ...
    'reserve_fraction',.3,'cumulative_charged_seconds',charged+seconds, ...
    'solver_sha256',solver_sha256,'new_boundary_solves',1, ...
    'purpose','cost-only coarse pilot' ,'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
if ~startsWith(id,'pilot'), receipt.purpose='current-version reference requalification and self-quadrature control'; end
save(file,'boundary','src','cfg_numeric','receipt','-v7.3');
closeout_s2b_json(fullfile(dest,['bie_' id '.json']),receipt); disp(receipt);
end
function ep=surface_derivative(x,A)
[~,ep]=closeout_s2b_surface(x,A);
end
