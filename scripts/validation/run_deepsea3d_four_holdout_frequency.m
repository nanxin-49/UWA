function receipt=run_deepsea3d_four_holdout_frequency(frequency_hz)
%RUN_DEEPSEA3D_FOUR_HOLDOUT_FREQUENCY Resume one locked frequency with a time cap.
arguments
    frequency_hz (1,1) double {mustBeMember(frequency_hz,[4500 5500 6500 7500])}
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_four_holdout_20261026');
if ~isfolder(out),mkdir(out);end
if ~isfile(fullfile(out,'four_holdout.mat'))
    validate_deepsea3d_four_holdout_frequency("init");
end
t=tic;
while true
    state=validate_deepsea3d_four_holdout_frequency("status");
    c=state.cases([state.cases.frequency_hz]==frequency_hz);
    if c.complete
        if ~(c.reference_ready && isfield(c.fk_guards,'gauss') && ...
                c.fk_guards.gauss>=1e-3 && isempty(fieldnames(c.fk_refinement)))
            break
        end
    end
    if toc(t)>67
        break
    end
    if c.complete
        stage="quadrature";
    elseif c.phase=="evaluation"
        stage="evaluate";
    elseif c.phase=="observe"
        stage="observe";
    elseif c.phase=="fk"
        stage="fk";
    else
        stage="next";
    end
    validate_deepsea3d_four_holdout_frequency(stage,frequency_hz);
end
state=validate_deepsea3d_four_holdout_frequency("status");
c=state.cases([state.cases.frequency_hz]==frequency_hz);
receipt=struct('frequency_hz',frequency_hz,'phase',c.phase,'complete',c.complete, ...
    'reference_ready',c.reference_ready,'G1_status',c.G1_status,'reason',c.reason, ...
    'elapsed_in_process_s',toc(t),'FE_runs',numel(c.runs),'FK_runs',numel(c.fk));
save(fullfile(out,sprintf('receipt_%d.mat',frequency_hz)),'receipt','-v7');
disp(receipt);
end
