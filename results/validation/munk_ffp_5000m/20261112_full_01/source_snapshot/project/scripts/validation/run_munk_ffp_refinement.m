function receipt=run_munk_ffp_refinement(stage,run_id)
%RUN_MUNK_FFP_REFINEMENT One targeted finer quadrature after a failed gate.
% Same 900s authorization, <=90s per execution; no physical-parameter change.
arguments
    stage (1,:) char
    run_id (1,:) char = '20261112_full_01'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'scripts','validation','support'),'-begin');
assert(isempty(regexp(run_id,'[^A-Za-z0-9_-]','once')));
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id);
auth=jsondecode(fileread(fullfile(d,'authorization.json'))); assert(auth.approved_total_budget_s>=900 && strcmp(auth.approved_by,'user'));
plan=load(fullfile(d,'frozen_plan.mat'),'plan'); cfg=plan.plan(6).config;
if strcmp(stage,'init')
    assert(~isfile(fullfile(d,'refinement_plan.mat')));
    cfg.panel_width_radpm=.0625; cfg.spectral_block=256;
    benchmark=jsondecode(fileread(fullfile(d,'refinement_matched_benchmark.json')));
    assert(benchmark.projected_12816_node_s+benchmark.extra_io_baseline_reserve_s<=90);
    reason='J03->J04 upper-column quadrature failed; same physics/finest layers, halve .125 to .0625. Block256 is separately checked for equivalence.';
    save(fullfile(d,'refinement_plan.mat'),'cfg','reason','benchmark','-v7');
    receipt=struct('stage','refinement_init','status','TARGETED_QUADRATURE_PREPARED','elapsed_s',0,'reason',reason,'predicted_wall_s',benchmark.projected_12816_node_s+6);
    local_json(fullfile(d,'refinement_init_receipt.json'),receipt); fprintf('%s\n',jsonencode(receipt)); return;
end
s=load(fullfile(d,'refinement_plan.mat'),'cfg'); proposed=s.cfg;
assert(any(strcmp(stage,{'block_control','solve'})));
if strcmp(stage,'solve'), tag='J07'; else, tag='block_control'; end
rf=fullfile(d,[tag '_receipt.json']); assert(~isfile(rf),'Preserve prior attempts; inspect before retry');
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=2*2^30,'munk_ffp:RefineRAM','Refinement block256 requires 2GiB starting RAM');
diary(fullfile(d,[tag '_diary.log'])); t=tic;
receipt=struct('stage',stage,'job_id',tag,'status','RUNNING','exit_status',NaN,'elapsed_s',0,'failed_or_unrun_checks',{{}});
try
    if strcmp(stage,'block_control')
        old=load(fullfile(d,'J06.mat'),'output');
        zr=[0 1000 3000 4900 4949.9 4950 4950.1 4999 5000]; xr=[-5 -.1 0 .1 5];
        [~,iz]=min(abs(old.output.depth_m-zr),[],1); [~,ix]=min(abs(old.output.range_m-xr.'),[],2);
        cfg.depth_m=old.output.depth_m(iz); cfg.range_m=old.output.range_m(ix); cfg.spectral_block=256;
        output=solve_layered_fluid_ffp(cfg); truth=old.output.green(iz,ix);
        valid=isfinite(truth); rel=norm(output.green(valid)-truth(valid))/norm(truth(valid));
        assert(rel<=1e-12,'Block-size regression fails: %.5g',rel);
        save(fullfile(d,'block_control.mat'),'output','truth','rel','-v7'); receipt.metrics=struct('block256_vs64_rel',rel);
    else
        control=jsondecode(fileread(fullfile(d,'block_control_receipt.json'))); assert(control.exit_status==0);
        output=solve_layered_fluid_ffp(proposed); save(fullfile(d,'J07.mat'),'output','-v7');
        receipt.metrics=struct('spectral_nodes',output.meta.spectral_node_count,'solve_wall_s',output.elapsed_s,'grid_size',size(output.green));
    end
    receipt.status='SAVED_NOT_YET_QUALIFIED'; receipt.exit_status=0;
catch err
    receipt.status='FAILED'; receipt.exit_status=1; receipt.failed_or_unrun_checks={err.identifier,err.message};
    fid=fopen(fullfile(d,[tag '_error.log']),'w'); fprintf(fid,'%s',getReport(err,'extended','hyperlinks','off')); fclose(fid);
end
receipt.elapsed_s=toc(t); local_json(rf,receipt); diary off; fprintf('%s\n',jsonencode(receipt));
end
function local_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); clean=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
