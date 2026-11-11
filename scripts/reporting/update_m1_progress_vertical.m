function output=update_m1_progress_vertical(run_dir)
%UPDATE_M1_PROGRESS_VERTICAL Merge this M1 stop without replacing other work.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
source=fullfile(root,'reports','project_progress.json'); before=sha(source);
receipt_path=fullfile(run_dir,'progress_update_receipt.json'); assert(~isfile(receipt_path));
d=jsondecode(fileread(source));
fid=fopen(fullfile(run_dir,'progress_source_before_final.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(d,PrettyPrint=true)); fclose(fid);
date=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd'));
d.updated_at=date; d.checked_at=date;
d.review_revision=[d.review_revision '; M1 resumed in attached R2025b Update 6/PID28836, then P8 RAM-reserve stop'];
d.headline=[d.headline ' 同日M1恢复后完成P4/P6及10组AS，但P8触发RAM门，仍为M1_INCOMPLETE。'];
d.overview{end+1}='同日后续M1：既有MCP版本/PID核实后，完成P4/P6和10组独立AS；已完成对照独立审计0失败。P8在RAM跌破启动可用量30%时中止，Z4/Z8和全部吸收层对照未派发，无生产推荐，不进入M2。';
d.review_scope=[d.review_scope ' 本次M1执行另核真实版本/PID、生产代码与相位v2，补默认关闭的每步验证观察；完成既有数据MAT/IEEE-754 CSV/SHA和参考/生产对照审计。未接管活塞计划或将其预算用于M1。'];
d.runtime_note=[d.runtime_note ' 同日后续M1使用现有本机MCP R2025b Update 6、PID28836、CPU double单线程串行。P4/P6与10组AS完成；P8用时552.6459804s，在入射976/1067步因可用RAM713232384B低于保留769405747B而停止。全部失败、先导重跑和审计计入14400s账；未使用SSH/GPU/新MATLAB。本次看板正式刷新结果以最新render_receipt为准，未做浏览器目视验收。'];
i=find(strcmp({d.workstreams.id},'broadband-reference')); assert(isscalar(i));
s=d.workstreams(i); s.checked_at=date; s.status='stopped'; s.status_label='M1_INCOMPLETE · P8 RAM门停止';
s.result_date=[date '（M1部分对照/资源停止）；正式S5/S6仍暂停'];
s.summary='M0封存保持。M1已完成P4/P6和10组AS，限定数值对照通过；P8因RAM门中止，步长与吸收层资格未完成，无生产推荐配置。';
s.done=[s.done ' 同日MCP恢复后执行生产P4/P6与10组AS，参考自身边缘、三频空间对照与4/6kHz PE-AS通过；独立收尾81行比较与MAT/CSV/SHA审计0失败。默认关闭的观察接口逐步记录三段传播，64²开关输出差0，已有M0宽带H消费者回放BER/SER=0。'];
s.meaning=[s.meaning ' 本轮证明指定均匀/平面Gaussian离散锚点的已完成对照，不形成完整M1配置资格。8kHz AS通过不能代替未完成生产PE或步长实测。'];
s.gap=[s.gap ' 当前P8仅直达1036/1036、入射977/1068、反射0/33平面；Z4/Z8及SP4/SP6/SP8均未派发。原吸收层既未判合格也未判不合格。'];
s.next='当前因RAM门收尾。后续用户继续指令后，先核现有会话/真实PID、输入与SHA、当前RAM和剩余预算，再设计不覆盖旧run的恢复；不自动重跑，不降低30%余量门，不进入M2。';
s.rerun=[s.rerun ' 当前MCP版本和PID已核实，原连接失败成为历史；本轮停止原因更新为P8 RAM门。原run02空diary与300s传输超时保留，经确认终态后以显式stdout捕获在新run03重跑P4且三表示圆盘差0。'];
new_refs={ ...
    'M1数值执行收尾与RAM停止','reports/m1_static_convergence_qualification_2026-11-11.md','report'; ...
    '已完成对照独立审计','results/validation/m1_convergence/20261111_m1_execution_03/partial_closeout_audit/checkpoint.json','audit'; ...
    '18项矩阵实际状态','results/validation/m1_convergence/20261111_m1_execution_03/matrix_status.csv','result'};
for j=1:size(new_refs,1)
    s.evidence(end+1)=struct('label',new_refs{j,1},'path',new_refs{j,2},'kind',new_refs{j,3},'revision','','expected_exists',true);
end
d.workstreams(i)=s;
event_id='m1-resumed-p8-ram-stop-20261111'; assert(~ismember(event_id,{d.events.id}));
d.events(end+1)=struct('id',event_id,'date',date,'title','M1恢复运行后在P8触发RAM预留门', ...
    'text','现有MCP确认R2025b Update 6/PID28836，完成P4/P6与10组AS；三频AS空间和4/6kHz生产同格对照及MAT/CSV/SHA独立审计通过。首次P4空diary/300s超时保留并计时，新run显式捕获stdout重跑且数值差0；首版参考审计错误保留后修正，仅重审已有MAT。P8在入射976步因RAM713232384B低于769405747B停止，Z4/Z8与全部SP未派发。故M1_INCOMPLETE，无生产推荐、不启动M2；下一步仅在继续指令和资源复核后恢复未完成M1。活塞计划及其他工作流记录/预算保留。', ...
    'workstream_id','broadband-reference');
d.next_actions(end+1)=struct('title','窗口M1停在RAM门，保留已完成证据', ...
    'text','P4/P6及全部AS已通过限定审计，P8部分轨迹保留；步长与吸收层尚未取得资格。', ...
    'completion','后续取得完整P8、Z4/Z8及阶段门后的吸收层证据，才可判M1最终资格。', ...
    'gate','当前不自动重跑或进入M2；保留30%启动可用RAM余量、14400s累计与3600s单作业限制。该预算不转用活塞计划。');
assert(strcmpi(before,sha(source)),'Progress source changed concurrently; merge again rather than overwrite.');
fid=fopen(source,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(d,PrettyPrint=true)); fclose(fid);
output=struct('status','M1_STOP_MERGED_OTHER_WORK_PRESERVED','source_before_sha256',before, ...
    'source_after_sha256',sha(source),'workstream','broadband-reference','event',event_id, ...
    'other_workstreams_preserved',numel(d.workstreams)-1,'existing_next_actions_preserved',numel(d.next_actions)-1);
fid=fopen(receipt_path,'w'); fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); fclose(fid);
fprintf('M1_PROGRESS %s\n',jsonencode(output));
end
function value=sha(file)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(file,'rb'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8'); if ~isempty(bytes), md.update(typecast(bytes,'int8')); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
