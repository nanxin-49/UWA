function output=validate_deepsea3d_preconditioner_decision()
%VALIDATE_DEEPSEA3D_PRECONDITIONER_DECISION Small diagnosis and strict RAM stop.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','deepsea_two_track_20261014');
file=fullfile(outdir,'preconditioner_decision.mat');
assert(~isfile(file),'Decision:Exists','Preserve prior checkpoint.');
timer=tic;[before,sys0]=memory;
cfg=struct('frequency_hz',600,'c0_mps',1500,'tx_m',[0 150.5 1500], ...
    'rx_m',[0 0 5],'center_m',[0 0],'profile',struct('height_m',0, ...
    'plateau_m',100,'support_m',101,'amplitudes_m',.01, ...
    'wavevectors_radpm',[pi/4 pi/4],'phases_rad',.2));
opts=struct('domain_halfwidth_m',4,'spacing_m',1,'vertical_degree',12,'top_depth_m',.8);
p=deepsea3d_fg_plan(cfg,opts);shape=[p.n^2,numel(p.s)];
action=@(x)reshape(deepsea3d_mapped_action(reshape(x,shape),p),[],1);
pre=@(x)deepsea3d_right_neumann_action(x,action);
right=@(x)action(pre(x));
v=exp(1i*pi/4*(p.xx+p.yy));truth=reshape(v(:)*sin(p.s.'+1),[],1);
other=conj(truth);rhs=action(truth);
linearity=norm(action(truth+2i*other)-action(truth)-2i*action(other))/norm(rhs);
polynomial=norm(right(truth)-(2*action(truth)-action(action(truth))))/norm(rhs);
[x0,flag0,res0,~,hist0]=gmres(action,rhs,8,1e-11,8);
[y,flag,res,~,hist]=gmres(right,rhs,8,1e-11,8);x=pre(y);
fresh=norm(action(x)-rhs)/norm(rhs);recover=norm(x-truth)/norm(truth);
% Explicit small matrix is an independent assembly of the same B columns.
d=numel(truth);matrix=complex(zeros(d));
for j=1:d
    e=zeros(d,1);e(j)=1;matrix(:,j)=action(e);
end
matrix_match=norm(matrix*truth-rhs)/norm(rhs);
explicit_recover=norm(matrix\rhs-truth)/norm(truth);
checks=struct('linearity',linearity,'polynomial',polynomial, ...
    'fresh_B_residual',fresh,'recover',recover,'matrix_match',matrix_match, ...
    'explicit_recover',explicit_recover,'baseline_recover',norm(x0-truth)/norm(truth));
assert(flag0==0 && flag==0 && max(struct2array_local(checks))<1e-8);
% Reuse the frozen raw/preconditioned residual regression without writing old artifacts.
[after,sys]=memory;estimate=663.19+4*18.5625;
eligible=estimate<.7*sys.PhysicalMemory.Available/2^20;
status="STOP_MEMORY_PREFLIGHT";
if eligible,status="RESOURCE_ELIGIBLE_NOT_RUN";end
baseline=load(fullfile(root,'results','validation','deepsea3d_mapped_direct','pm_restart12.mat'),'output');
output=struct('status',status,'route','right P=2I-B; BP=2B-B^2', ...
    'small_checks',checks,'gmres_flag',flag,'small_preconditioned_relres',res, ...
    'baseline_small_relres',res0,'history',hist,'baseline_small_history',hist0, ...
    'available_before_MiB',sys0.PhysicalMemory.Available/2^20, ...
    'available_after_MiB',sys.PhysicalMemory.Available/2^20, ...
    'required_estimate_MiB',estimate,'resource_ceiling_MiB',.7*sys.PhysicalMemory.Available/2^20, ...
    'process_before_MiB',before.MemUsedMATLAB/2^20,'process_after_MiB',after.MemUsedMATLAB/2^20, ...
    'peak_memory_MiB',NaN,'peak_kind','unknown; snapshots only', ...
    'new_PM_solve_executed',false,'linear_solved_PM',false, ...
    'reference_candidate',false,'reference_ready',false, ...
    'G1_target_PM','OPEN','WP2_formal',false,'wall_s',toc(timer), ...
    'baseline',baseline.output,'manufactured_cfg',cfg,'manufactured_options',opts);
save(file,'output','-v7');
summary=rmfield(output,'baseline');
summary.baseline_residuals=struct('fresh_B',baseline.output.fresh_preconditioned_residual, ...
    'discrete',baseline.output.discrete,'independent4',baseline.output.independent4);
write_json(fullfile(outdir,'preconditioner_decision.json'),summary);
fprintf('A %s; checks %.3g; required %.2f / ceiling %.2f MiB; wall %.3fs\n', ...
    status,max(struct2array_local(checks)),estimate,output.resource_ceiling_MiB,output.wall_s);
end
function values=struct2array_local(s)
values=cell2mat(struct2cell(s));
end
function write_json(file,data)
fid=fopen(file,'w');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(data,PrettyPrint=true));
end
