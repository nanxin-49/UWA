function output = render_project_progress_dashboard(overrides)
%RENDER_PROJECT_PROGRESS_DASHBOARD Render the maintained progress record.
%   setup_vertical_project; render_project_progress_dashboard();
%   Optional fields: source_path, output_dir, markdown_path (absolute paths).
%   Only JSON/text and exact referenced-file existence are read. No acoustic
%   solver, MAT/FIG loader, Git command, network, or numerical audit is run.
%   HTML/Markdown are replaceable views; reports/project_progress.json is the
%   sole maintained record. Existing views are intentionally refreshed.
arguments
    overrides (1,1) struct = struct()
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg = struct('source_path',fullfile(root,'reports','project_progress.json'), ...
    'output_dir',fullfile(root,'results','visualization','project_progress'), ...
    'markdown_path',fullfile(root,'reports','project_progress.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    assert(isfield(cfg,names{ii}),'progress:UnknownOption', ...
        'Unknown option: %s',names{ii});
    cfg.(names{ii}) = char(overrides.(names{ii}));
end
local_bound_path(root,cfg.source_path,'reports');
local_bound_path(root,cfg.markdown_path,'reports');
local_bound_path(root,cfg.output_dir,'results/visualization');
assert(endsWith(cfg.source_path,'.json') && ...
    endsWith(cfg.markdown_path,'.md'),'progress:Extension', ...
    'Expected a JSON source and a Markdown output.');
template_path = fullfile(root,'scripts','reporting','templates', ...
    'project_progress.html');
record = jsondecode(fileread(cfg.source_path));
local_validate_record(record);

checks = struct('workstream_id',{},'path',{},'exists',{}, ...
    'expected_exists',{},'kind',{},'revision',{});
for ii = 1:numel(record.workstreams)
    stream = record.workstreams(ii);
    for jj = 1:numel(stream.evidence)
        ref = stream.evidence(jj);
        local_validate_relative(ref.path);
        checks(end+1) = struct('workstream_id',stream.id, ...
            'path',ref.path,'exists',isfile(fullfile(root,ref.path)), ...
            'expected_exists',logical(ref.expected_exists), ...
            'kind',ref.kind,'revision',ref.revision); %#ok<AGROW>
    end
end
generated_at = char(datetime('now','TimeZone','Asia/Shanghai', ...
    'Format','yyyy-MM-dd HH:mm:ss Z'));
root_href = local_relative_path(cfg.output_dir,root);
markdown_root_href = local_relative_path(fileparts(cfg.markdown_path),root);
payload = struct('record',record,'file_checks',checks, ...
    'root_href',root_href,'markdown_root_href',markdown_root_href, ...
    'generated_at',generated_at);
% An embedded JSON script must never be terminated by source prose.
payload_json = strrep(jsonencode(payload),'<','\u003c');
html = fileread(template_path);
html = strrep(html,'@@PROGRESS_JSON@@',payload_json);
html = strrep(html,'@@ROOT_HREF@@',root_href);
markdown_href = [root_href local_escape_path( ...
    local_relative_path(root,fileparts(cfg.markdown_path))) ...
    local_escape_path(local_filename(cfg.markdown_path))];
html = strrep(html,'@@MARKDOWN_HREF@@',markdown_href);
assert(~contains(html,'@@'),'progress:Template', ...
    'Unexpanded template token.');
markdown = local_markdown(record,checks,markdown_root_href);
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
html_path = fullfile(cfg.output_dir,'index.html');
local_write_utf8(html_path,html);
local_write_utf8(cfg.markdown_path,markdown);
unexpected = checks([checks.exists] ~= [checks.expected_exists]);
status = 'GENERATED_FILE_CHECKS_PASS';
if ~isempty(unexpected), status = 'GENERATED_WITH_REFERENCE_WARNINGS'; end
receipt = struct('status',status,'generated_at',generated_at, ...
    'matlab_release',version,'source_updated_at',record.updated_at, ...
    'review_checked_at',record.checked_at,'file_checks',checks, ...
    'failed_checks',numel(unexpected),'unexpected_references',unexpected, ...
    'source_sha256',local_sha256(cfg.source_path), ...
    'renderer_sha256',local_sha256([mfilename('fullpath') '.m']), ...
    'template_sha256',local_sha256(template_path), ...
    'html_sha256',local_sha256(html_path), ...
    'markdown_sha256',local_sha256(cfg.markdown_path), ...
    'new_acoustic_computations',0, ...
    'audit_scope','Exact-file existence only; historical numerical audits not rerun');
receipt_path = fullfile(cfg.output_dir,'render_receipt.json');
local_write_utf8(receipt_path,jsonencode(receipt,PrettyPrint=true));
output = struct('status',status,'html_path',html_path, ...
    'markdown_path',cfg.markdown_path,'receipt_path',receipt_path, ...
    'checked_references',numel(checks),'failed_checks',numel(unexpected));
fprintf('PROGRESS_RENDER status=%s refs=%d unexpected=%d\n', ...
    status,numel(checks),numel(unexpected));
fprintf('HTML=%s\nMARKDOWN=%s\nRECEIPT=%s\n', ...
    html_path,cfg.markdown_path,receipt_path);
end

function local_validate_record(d)
required = {'schema_version','title','subtitle','updated_at','checked_at', ...
    'review_revision','headline','overview','reading_note','review_scope', ...
    'runtime_note','workstreams','events','next_actions','glossary','maintenance'};
assert(all(isfield(d,required)) && d.schema_version == 1, ...
    'progress:Schema','Incomplete record or unsupported schema_version.');
ids = {d.workstreams.id};
assert(numel(unique(ids)) == numel(ids),'progress:DuplicateId', ...
    'Workstream IDs must be unique.');
allowed = {'closed','partial','paused','pending','stopped'};
fields = {'name','status','status_label','result_date','checked_at','summary','why', ...
    'done','meaning','gap','next','details','saved_results', ...
    'historical_audit','rerun','evidence'};
for ii = 1:numel(d.workstreams)
    s = d.workstreams(ii);
    assert(all(isfield(s,fields)) && ismember(s.status,allowed), ...
        'progress:Workstream','Incomplete workstream or invalid status.');
    assert(~isempty(regexp(s.id,'^[a-z][a-z0-9-]*$','once')), ...
        'progress:Id','Workstream IDs must be safe HTML anchors.');
    assert(all(isfield(s.evidence, ...
        {'label','path','kind','revision','expected_exists'})), ...
        'progress:Evidence','Incomplete evidence reference.');
end
assert(numel(unique({d.events.id})) == numel(d.events), ...
    'progress:EventId','Event IDs must be unique.');
for ii = 1:numel(d.events)
    assert(ismember(d.events(ii).workstream_id,ids), ...
        'progress:EventTarget','Event refers to an unknown workstream.');
end
local_validate_relative(d.maintenance.source);
end

function local_validate_relative(p)
p = strrep(char(p),'\','/');
parts = strsplit(p,'/');
assert(~isempty(p) && ~startsWith(p,'/') && ~contains(p,':') && ...
    ~any(ismember(parts,{'..','.',''})) && ...
    ~strcmpi(parts{1},'cash'),'progress:UnsafeReference', ...
    'Evidence must be a project-relative path outside the retired archive.');
end

function local_bound_path(root,p,category)
p = strrep(char(p),'\','/');
base = [strrep(root,'\','/') '/' category '/'];
assert(startsWith(lower([p '/']),lower(base)) && ...
    ~any(ismember(strsplit(p,'/'),{'..','.'})), ...
    'progress:OutputBoundary','Path must stay inside %s.',category);
end

function href = local_relative_path(from_dir,to_dir)
from = strsplit(strrep(from_dir,'\','/'),'/');
to = strsplit(strrep(to_dir,'\','/'),'/');
common = 0;
while common < min(numel(from),numel(to)) && ...
        strcmpi(from{common+1},to{common+1})
    common = common+1;
end
assert(common > 0,'progress:PathRoot','Paths must share a filesystem root.');
href = [repmat('../',1,numel(from)-common) ...
    strjoin(to(common+1:end),'/')];
if ~isempty(href) && ~endsWith(href,'/'), href = [href '/']; end
end

function name = local_filename(path)
[~,stem,ext] = fileparts(path);
name = [stem ext];
end

function text = local_escape_path(text)
text = strrep(text,'%','%25');
text = strrep(text,' ','%20');
text = strrep(text,'#','%23');
text = strrep(text,'?','%3F');
text = strrep(text,'&','%26');
text = strrep(text,'"','%22');
text = strrep(text,'''','%27');
text = strrep(text,'<','%3C');
text = strrep(text,'>','%3E');
text = strrep(text,'(','%28');
text = strrep(text,')','%29');
end

function md = local_markdown(d,checks,root_href)
lines = [ ...
    "<!-- GENERATED: edit reports/project_progress.json, then run render_project_progress_dashboard. -->"; ...
    "# " + string(d.title); ""; string(d.subtitle); ""; ...
    "记录更新：" + string(d.updated_at) + "；本次核查：" + string(d.checked_at) + "。"; ...
    ""; string(d.headline); ""; "## 现在在哪里"; ""];
for ii = 1:numel(d.overview)
    lines = [lines; string(d.overview{ii}); ""]; %#ok<AGROW>
end
lines = [lines; string(d.reading_note); ""; "## 各条工作线"; ""];
keys = {'why','done','meaning','gap','next','details', ...
    'saved_results','historical_audit','rerun'};
labels = {'为什么做','做了什么','说明什么','还差什么','下一步', ...
    '技术范围','已有结果','历史审计','当前重跑条件'};
for ii = 1:numel(d.workstreams)
    s = d.workstreams(ii);
    lines = [lines; "### " + string(s.name); ""; ...
        "状态：" + string(s.status_label) + "；结果日期：" + ...
        string(s.result_date) + "；本次核查：" + string(s.checked_at) + "。"; ...
        ""; string(s.summary); ""]; %#ok<AGROW>
    for jj = 1:numel(keys)
        lines = [lines; "**" + string(labels{jj}) + "：**" + ...
            string(s.(keys{jj})); ""]; %#ok<AGROW>
    end
    lines = [lines; "证据入口："; ""]; %#ok<AGROW>
    for jj = 1:numel(s.evidence)
        e = s.evidence(jj);
        k = find(strcmp({checks.workstream_id},s.id) & ...
            strcmp({checks.path},e.path),1);
        if checks(k).exists
            line = "- [" + string(e.label) + "](" + string(root_href) + ...
                string(local_escape_path(e.path)) + ")";
        else
            line = "- **" + string(e.label) + "：当前缺失**";
        end
        line = line + " · `" + string(e.path) + "`";
        if ~isempty(e.revision)
            line = line + " · 历史版本 `" + string(e.revision) + "`";
        end
        if checks(k).exists ~= e.expected_exists
            line = line + " · 与记录预期不符，待复核";
        end
        lines = [lines; line]; %#ok<AGROW>
    end
    lines = [lines; ""]; %#ok<AGROW>
end
lines = [lines; "## 最近发生了什么"; ""];
for ii = 1:numel(d.events)
    e = d.events(ii);
    lines = [lines; "### " + string(e.date) + " · " + ...
        string(e.title); ""; string(e.text); ""]; %#ok<AGROW>
end
lines = [lines; "## 重新开始先做什么"; ""];
for ii = 1:numel(d.next_actions)
    a = d.next_actions(ii);
    lines = [lines; "### " + string(ii) + ". " + string(a.title); ...
        ""; string(a.text); ""; "完成标志：" + string(a.completion); ...
        ""; "执行边界：" + string(a.gate); ""]; %#ok<AGROW>
end
lines = [lines; "## 术语随手查"; ""];
for ii = 1:numel(d.glossary)
    g = d.glossary(ii);
    lines = [lines; "**" + string(g.term) + "：**" + ...
        string(g.meaning); ""]; %#ok<AGROW>
end
lines = [lines; "## 记录与核查范围"; ""; string(d.review_scope); ...
    ""; string(d.runtime_note); ""; ...
    "本次核查代码基线：`" + string(d.review_revision) + "`。"; ""; ...
    "唯一维护源：`" + string(d.maintenance.source) + "`。"; ""; ...
    string(d.maintenance.rule); ""; string(d.maintenance.qualification); ...
    ""; "```matlab"; string(d.maintenance.command); "```"; ""];
md = char(join(lines,newline));
end

function local_write_utf8(path,text)
fid = fopen(path,'w','n','UTF-8');
assert(fid > 0,'progress:Write','Cannot write %s.',path);
guard = onCleanup(@()fclose(fid));
fprintf(fid,'%s',text);
end

function value = local_sha256(path)
fid = fopen(path,'rb');
assert(fid > 0,'progress:Read','Cannot hash %s.',path);
guard = onCleanup(@()fclose(fid));
bytes = fread(fid,Inf,'*uint8');
% Only the bundled Java SHA-256 primitive; no custom classes or classpaths.
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes);
digest = typecast(engine.digest(),'uint8');
value = lower(reshape(dec2hex(digest,2).',1,[]));
end
