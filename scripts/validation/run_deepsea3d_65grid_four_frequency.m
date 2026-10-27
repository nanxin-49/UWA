function receipt=run_deepsea3d_65grid_four_frequency(frequency_hz)
%RUN_DEEPSEA3D_65GRID_FOUR_FREQUENCY One prelocked frequency; no retry/refinement.
arguments
    frequency_hz (1,1) double {mustBeMember(frequency_hz,[4750 7250 4125 6562.5])}
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_65grid_four_20261027');
if ~isfolder(out),mkdir(out);end
if ~isfile(fullfile(out,'65grid_four.mat'))
    validate_deepsea3d_65grid_four_frequency("init");
end
t=tic;
while true
    state=validate_deepsea3d_65grid_four_frequency("status");
    c=state.cases([state.cases.frequency_hz]==frequency_hz);
    if c.complete
        break
    end
    if toc(t)>120
        break
    end
    if c.phase=="evaluation"
        stage="evaluate";
    elseif c.phase=="observe"
        stage="observe";
    elseif c.phase=="fk"
        stage="fk";
    else
        stage="next";
    end
    validate_deepsea3d_65grid_four_frequency(stage,frequency_hz);
end
state=validate_deepsea3d_65grid_four_frequency("status");
c=state.cases([state.cases.frequency_hz]==frequency_hz);
receipt=struct('frequency_hz',frequency_hz,'phase',c.phase,'complete',c.complete, ...
    'reference_ready',c.reference_ready,'G1_status',c.G1_status,'reason',c.reason, ...
    'elapsed_in_process_s',toc(t),'FE_runs',numel(c.runs),'FK_runs',numel(c.fk));
save(fullfile(out,sprintf('receipt_%g.mat',frequency_hz)),'receipt','-v7');
disp(receipt);
end
