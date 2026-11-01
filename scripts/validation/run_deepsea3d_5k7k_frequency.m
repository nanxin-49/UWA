function receipt=run_deepsea3d_5k7k_frequency(frequency_hz)
%RUN_DEEPSEA3D_5K7K_FREQUENCY Resume one locked frequency with a time cap.
arguments
    frequency_hz (1,1) double {mustBeMember(frequency_hz,[5000 7000])}
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_5k7k_20261026');
if ~isfolder(out),mkdir(out);end
if ~isfile(fullfile(out,'five_seven.mat'))
    validate_deepsea3d_5k7k_frequency("init");
end
t=tic;
while true
    state=validate_deepsea3d_5k7k_frequency("status");
    c=state.cases([state.cases.frequency_hz]==frequency_hz);
    if c.complete && ~(c.reference_ready && c.fk_guards.gauss>=1e-3 && isempty(fieldnames(c.fk_refinement)))
        break
    end
    if toc(t)>235
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
    validate_deepsea3d_5k7k_frequency(stage,frequency_hz);
end
state=validate_deepsea3d_5k7k_frequency("status");
c=state.cases([state.cases.frequency_hz]==frequency_hz);
receipt=struct('frequency_hz',frequency_hz,'phase',c.phase,'complete',c.complete, ...
    'reference_ready',c.reference_ready,'G1_status',c.G1_status,'reason',c.reason, ...
    'elapsed_in_process_s',toc(t),'FE_runs',numel(c.runs),'FK_runs',numel(c.fk));
save(fullfile(out,sprintf('receipt_%d.mat',frequency_hz)),'receipt','-v7');
disp(receipt);
end
