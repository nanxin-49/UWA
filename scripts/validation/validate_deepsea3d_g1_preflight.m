function results=validate_deepsea3d_g1_preflight()
%VALIDATE_DEEPSEA3D_G1_PREFLIGHT Small bounded reference feasibility probe.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_g1_finite');
if ~isfolder(out), mkdir(out); end
cfg=struct('frequency_hz',6000,'c0_mps',1500,'tx_m',[0 0 1500], ...
    'rx_m',[0 0 5],'center_m',[0 0],'half_width_m',16,'spacing_m',.025, ...
    'quadrature','trap','profile',struct('height_m',0,'plateau_m',8,'support_m',12, ...
    'amplitudes_m',[.02 .01],'wavevectors_radpm',[.4 0;.2 .6],'phases_rad',[.2 .7]));
rows=struct([]); results=struct('cfg',cfg,'runs',struct([]),'decision','UNASSESSED');
for spacing=[.25 .125]
    opts=struct('domain_halfwidth_m',24,'spacing_m',spacing,'order',6,'residual_order',20);
    a=deepsea3d_dirichlet_fe(cfg,opts);
    row=struct('spacing_m',a.spacing_m,'nodes',a.nodes,'wall_s',a.wall_s,'cpu_s',a.cpu_s, ...
        'correction_real',real(a.rough_correction),'correction_imag',imag(a.rough_correction), ...
        'boundary_residual',a.boundary_residual_relative, ...
        'order_4_6_change',abs(a.partials(6)-a.partials(4))/abs(a.partials(6)), ...
        'stored_spectral_MiB',a.stored_spectral_MiB,'max_eta_kz',a.max_eta_kz);
    if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
    if isempty(results.runs), results.runs=a; else, results.runs(end+1)=a; end
    save(fullfile(out,'preflight.mat'),'results','-v7');
end
results.rows=rows;
results.decision='READ_PRECHECK_BEFORE_REFERENCE_CAMPAIGN';
save(fullfile(out,'preflight.mat'),'results','-v7');
writetable(struct2table(rows),fullfile(out,'preflight.csv')); disp(struct2table(rows));
end
