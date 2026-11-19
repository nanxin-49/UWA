function k1_recovery_event(root,label,state)
%K1_RECOVERY_EVENT Durable per-step observability, no physical arithmetic.
event=struct('label',label,'state',state,'PID',feature('getpid'), ...
    'time',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SSS Z')));
p=fullfile(root,'recovery_events.jsonl'); fid=fopen(p,'a'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(event)); fclose(fid);
fprintf('K1_RECOVERY_EVENT %s %s\n',label,state);
control=fullfile(fileparts(root),'control'); if isfolder(control)
    fid=fopen(fullfile(control,'model_event.txt'),'w'); assert(fid>=0); fprintf(fid,'%s %s\n',label,state); fclose(fid);
end
end
