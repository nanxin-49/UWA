function receipt=run_munk_ffp_preflight(run_id)
%RUN_MUNK_FFP_PREFLIGHT Existing baselines, small layered controls and cost pilot.
% No full-field solve, no Bellhop/PE, no broadband, no automatic budget increase.
arguments
    run_id (1,:) char = '20261112_preflight_01'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(isempty(regexp(run_id,'[^A-Za-z0-9_-]','once')));
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id);
assert(isfolder(d),'Save design/input identity first.');
assert(~isfile(fullfile(d,'preflight_receipt.json')),'Never overwrite an attempted preflight.');
addpath(fullfile(root,'scripts','validation','support'),'-begin'); diary(fullfile(d,'preflight_diary.log')); t=tic;
receipt=struct('run_id',run_id,'status','RUNNING','elapsed_s',0,'exit_status',NaN, ...
    'execution_mode','existing MATLAB MCP only','matlab_release',version,'pid',feature('getpid'), ...
    'full_field_solved',false,'failed_or_unrun_checks',{{}});
try
    cfg=make_munk_5000m_ffp_config(10); cfg.layer_edges_m=(0:2:20).';
    cfg.layer_c_mps=1500*ones(10,1); cfg.source_depth_m=12; cfg.water_depth_m=20;
    cfg.bottom_type='outgoing'; cfg.top_type='outgoing'; cfg.subtract_reference=false;
    cfg.range_m=[-5 -1 -.1 0 .1 1 5]; cfg.depth_m=[.25;2;5;10;14;19];
    cfg.max_wall_s=20; kv=[0 .2 5 39 45 80]-.2i;
    [g,~]=layered_fluid_green_kernel(cfg,kv,cfg.depth_m);
    q=sqrt(complex((2*pi*1e4/1500)^2-kv.^2));
    exact=(1i./(2*q)).*exp(1i*abs(cfg.depth_m-12).*q);
    kernel_error=norm(g-exact,'fro')/norm(exact,'fro');
    assert(kernel_error<=1e-10,'munk_ffp:UniformKernel','Uniform kernel error %.6g',kernel_error);
    delta=1e-5; zz=12+[-2 -1 0 1 2]*delta;
    [jfield,~]=layered_fluid_green_kernel(cfg,kv,zz(:));
    dp=(-3*jfield(3,:)+4*jfield(4,:)-jfield(5,:))/(2*delta);
    dm=(3*jfield(3,:)-4*jfield(2,:)+jfield(1,:))/(2*delta);
    jump_error=max(abs(dp-dm+1)); assert(jump_error<2e-6,'Source jump failed');
    uniform_free=solve_layered_fluid_ffp(cfg);
    [zz,xx]=ndgrid(cfg.depth_m,cfg.range_m); points=[xx(:),zeros(numel(xx),1),zz(:)];
    truth=reshape(deepsea3d_point_field(points,[0 0 12],2*pi*1e4/1500),size(xx));
    free_rel=norm(uniform_free.green-truth,'fro')/norm(truth,'fro');
    free_max_phase=max(abs(angle(uniform_free.green.*conj(truth))),[],'all');
    save(fullfile(d,'uniform_free.mat'),'uniform_free','truth','-v7');
    cfg.top_type='pressure_release'; uniform_release=solve_layered_fluid_ffp(cfg);
    mirror=complex(zeros(size(truth)));
    for j=1:numel(mirror)
        p=struct('frequency_hz',1e4,'c0_mps',1500,'tx_m',[0 0 12],'rx_m',points(j,:));
        old=deepsea3d_flat(p,0); mirror(j)=old.total_native;
    end
    release_rel=norm(uniform_release.green-mirror,'fro')/norm(mirror,'fro');
    save(fullfile(d,'uniform_release.mat'),'uniform_release','mirror','-v7');
    [~,weyl_meta]=weyl_point_source_reference_vertical([0 .1 1],2,1e4,1500);
    assert(free_rel<=1e-3 && free_max_phase<=.01 && release_rel<=1e-3,'munk_ffp:UniformIntegral', ...
        'Uniform transform fails: free %.4g phase %.4g mirror %.4g',free_rel,free_max_phase,release_rel);
    points_z=[0;1000;3000;4900;4949.9;4950;4950.1;4999;5000];
    p=make_munk_5000m_ffp_config(10); p.range_m=[-5 -.1 0 .1 5]; p.depth_m=points_z;
    p.max_wall_s=30; coarse=solve_layered_fluid_ffp(p); save(fullfile(d,'munk_pilot_coarse.mat'),'coarse','-v7');
    p.panel_width_radpm=.25; fine=solve_layered_fluid_ffp(p); save(fullfile(d,'munk_pilot_fine.mat'),'fine','-v7');
    mask=~fine.source_singularity_mask;
    integral_change=norm(coarse.green(mask)-fine.green(mask))/norm(fine.green(mask));
    axis=fine.range_m==0; valid=~fine.source_singularity_mask(:,axis);
    assert(all(isfinite(fine.green(valid,axis))));
    surface_residual=max(abs(fine.green(fine.depth_m==0,:)),[],'all')/max(abs(fine.green(mask)));
    % Explicit temporaries avoid unsupported indexing of an expression.
    mirrored=fliplr(fine.green); symmetry=norm(fine.green(mask)-mirrored(mask))/norm(fine.green(mask));
    tau=munk_vertical_travel_time_reference(points_z,4950,p.munk_parameters);
    time_rows=struct([]);
    for layer_h=[10 5 2.5]
        model=make_munk_5000m_ffp_config(layer_h); tl=zeros(size(points_z));
        for j=1:numel(points_z)
            lo=min(points_z(j),4950); hi=max(points_z(j),4950);
            lengths=max(0,min(model.layer_edges_m(2:end),hi)-max(model.layer_edges_m(1:end-1),lo));
            tl(j)=sum(lengths./model.layer_c_mps);
        end
        time_rows=[time_rows struct('layer_h_m',layer_h,'max_travel_error_s',max(abs(tl-tau)), ...
            'max_onepass_phase_error_rad',2*pi*1e4*max(abs(tl-tau)))]; %#ok<AGROW>
    end
    writetable(struct2table(time_rows),fullfile(d,'vertical_delay_layer_convergence.csv'));
    save(fullfile(d,'vertical_delay_reference.mat'),'tau','points_z','time_rows','-v7');
    % One actual spectral block on the requested combined receiver grid, not a field result.
    assert(toc(t)<120,'Preflight budget reserve');
    zbench=unique([(0:5:5000)';(4900:.01:5000)']); xbench=-5:.05:5;
    kbench=linspace(.1,4,64)-.2i; bt=tic; gb=layered_fluid_green_kernel(p,kbench,zbench);
    dummy=gb*besselj(0,kbench.'*abs(xbench)); block_s=toc(bt); clear gb dummy;
    pilot=struct('z_count',numel(zbench),'x_count',numel(xbench),'spectral_block',64, ...
        'measured_kernel_and_transform_block_s',block_s,'estimated_full_spectral_nodes',fine.meta.spectral_node_count, ...
        'linear_block_extrapolation_s',block_s*ceil(fine.meta.spectral_node_count/64), ...
        'cost_limitations','One block, no subtraction/save/render costs; use conservative factor and reserve.');
    fid=fopen(fullfile(d,'cost_pilot.json'),'w'); fprintf(fid,'%s',jsonencode(pilot,'PrettyPrint',true)); fclose(fid);
    receipt.metrics=struct('uniform_kernel_rel',kernel_error,'source_jump_max_abs',jump_error, ...
        'uniform_free_complex_rel',free_rel,'uniform_free_max_phase_rad',free_max_phase, ...
        'uniform_mirror_complex_rel',release_rel,'existing_weyl_rel',weyl_meta.relative_error_to_green, ...
        'munk_pilot_integral_refinement_rel',integral_change,'munk_surface_rel',surface_residual, ...
        'munk_symmetry_rel',symmetry,'cost_pilot',pilot);
    receipt.status='PREFLIGHT_COMPLETE_FULL_FIELD_PENDING'; receipt.exit_status=0;
catch err
    receipt.status='PREFLIGHT_FAILED'; receipt.exit_status=1;
    receipt.failed_or_unrun_checks={err.identifier,err.message};
    fid=fopen(fullfile(d,'error.log'),'w'); fprintf(fid,'%s',getReport(err,'extended','hyperlinks','off')); fclose(fid);
end
receipt.elapsed_s=toc(t); assert(receipt.elapsed_s<180,'Preflight exceeded budget');
fid=fopen(fullfile(d,'preflight_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,'PrettyPrint',true)); fclose(fid); diary off;
fprintf('%s\n',jsonencode(receipt));
if receipt.exit_status~=0, error('munk_ffp:PreflightFailed','Read saved error and receipt.'); end
end
