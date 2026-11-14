function receipt = validate_bellhop_rotated_cost_pilot(executable,output_dir)
%VALIDATE_BELLHOP_ROTATED_COST_PILOT Bounded cost/capacity, not field acceptance.
% Four components only, <=60s total and <=20s each; no complete 2D grid.
arguments
    executable (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
tags={'long_Q625','long_Q3125','near_q11','near_q201'};
rows=struct('tag',{},'points',{},'beams',{},'Q_spacing_m',{},'max_steps',{}, ...
    'wall_s',{},'process_s',{},'sampled_peak_bytes',{},'range_count',{},'q_count',{});
receipt=struct('status','FAILED','result_dir',output_dir,'max_wall_s',60,'max_job_s',20, ...
    'scope','Four sparse cost/capacity components only; no complete-grid field or convergence qualification');
try
    q_counts=[11 201];
    for j=1:4
        assert(toc(t)<60,'Cost pilot cumulative time gate.');
        if j<=2, z=0:5:5000; x=[-5 0 5]; else, z=4900:.25:5000; x=linspace(-5,5,q_counts(j-2)); end
        cfg=make_bellhop_rotated_multifield_config('up',16,z,x);
        cfg.ssp_spacing_m=.625; if j==2, cfg.ssp_spacing_m=.3125; end
        cfg.step_m=.5; cfg.beam_count=51;
        fan=1.1*rad2deg(atan(1.05*max(abs(x))/min(cfg.receiver_range_m)))+.001;
        cfg.angle_limits_deg=[-fan fan]; cfg.range_limit_m=max(cfg.receiver_range_m);
        a=run_bellhop_rotated_field_component(executable,cfg,fullfile(output_dir,tags{j}), ...
            max_seconds=min(20,60-toc(t)));
        rows(end+1)=struct('tag',tags{j},'points',numel(z)*numel(x),'beams',cfg.beam_count, ...
            'Q_spacing_m',cfg.ssp_spacing_m,'max_steps',max(a.wall(:,6)), ...
            'wall_s',a.receipt.elapsed_s,'process_s',a.receipt.process.elapsed_s, ...
            'sampled_peak_bytes',a.receipt.process.peak_working_set_bytes, ...
            'range_count',numel(z),'q_count',numel(x)); %#ok<AGROW>
        clear a
    end
    receipt.status='COMPLETE_COST_PILOT_ONLY';
catch e
    receipt.error=getReport(e,'extended','hyperlinks','off');
end
receipt.elapsed_s=toc(t); receipt.rows=rows;
if ~isempty(rows), writetable(struct2table(rows),fullfile(output_dir,'timings.csv')); end
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
save(fullfile(output_dir,'result.mat'),'receipt','rows'); disp(jsonencode(receipt));
end