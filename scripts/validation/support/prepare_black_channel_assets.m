function assets = prepare_black_channel_assets(project_root)
%PREPARE_BLACK_CHANNEL_ASSETS Fetch only BLACK data and pinned official code.
% Existing raw files must match the official digest; they are never replaced.
if nargin == 0
    project_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
case_dir = fullfile(project_root, 'results', 'validation', 'black_channel_basics');
raw_dir = fullfile(case_dir, 'raw');
revision = 'af1d3814af9199f1e420a248ab27edaced5339cf';
vendor_dir = fullfile(case_dir, 'vendor', revision);
if ~isfolder(raw_dir), mkdir(raw_dir); end
if ~isfolder(vendor_dir), mkdir(vendor_dir); end
manifest_file = fullfile(raw_dir, 'manifest.json');
if isfile(manifest_file)
    manifest = jsondecode(fileread(manifest_file));
    assert(strcmp(manifest.source_revision, revision), 'black:Revision', ...
        'Asset manifest has a different source revision.');
else
    record = webread('https://zenodo.org/api/records/21287414', weboptions('Timeout', 60));
    manifest = struct('record_id', record.id, 'doi', record.doi, ...
        'dataset_version', record.metadata.version, 'license', record.metadata.license.id, ...
        'source_revision', revision, 'created_utc', char(datetime('now','TimeZone','UTC')), ...
        'files', repmat(struct('name','','size',0,'official_checksum','','url','','sha256',''),1,2));
    names = {'black.mat', 'black_noise.mat'};
    for ii = 1:numel(names)
        index = find(strcmp({record.files.key}, names{ii}));
        assert(isscalar(index), 'black:MissingAsset', 'BLACK file missing from record.');
        item = record.files(index);
        manifest.files(ii) = struct('name', names{ii}, 'size', item.size, ...
            'official_checksum', item.checksum, 'url', item.links.self, 'sha256', '');
    end
end
for ii = 1:numel(manifest.files)
    item = manifest.files(ii);
    assert(any(strcmp(item.name, {'black.mat','black_noise.mat'})), 'black:Scope', 'BLACK files only.');
    filename = fullfile(raw_dir, item.name);
    if ~isfile(filename)
        part = [filename '.part'];
        assert(~isfile(part), 'black:PartialDownload', 'Inspect existing partial download: %s', part);
        fprintf('Downloading %s (%d bytes) ...\n', item.name, item.size);
        websave(part, item.url, weboptions('Timeout', 180));
        info = dir(part);
        assert(info.bytes == item.size, 'black:DownloadSize', 'Wrong download size.');
        expected = extractAfter(string(item.official_checksum), 'md5:');
        assert(strcmp(black_file_hash(part, 'MD5'), expected), 'black:DownloadDigest', 'MD5 mismatch.');
        movefile(part, filename);
    end
    info = dir(filename);
    assert(info.bytes == item.size, 'black:AssetSize', 'Existing raw size mismatch.');
    assert(strcmp(['md5:' black_file_hash(filename,'MD5')], item.official_checksum), ...
        'black:AssetDigest', 'Existing raw checksum mismatch: %s', item.name);
    manifest.files(ii).sha256 = black_file_hash(filename);
    fprintf('Verified %s MD5 and SHA-256\n', item.name);
end
code_names = {'replay.m','unpack.m','noisegen.m','LICENSE'};
previous_code = [];
if isfield(manifest,'code'), previous_code = manifest.code; end
manifest.code = repmat(struct('name','','url','','sha256',''),1,numel(code_names));
for ii = 1:numel(code_names)
    filename = fullfile(vendor_dir, code_names{ii});
    if strcmp(code_names{ii}, 'LICENSE'), relative = 'LICENSE';
    else, relative = ['src/' code_names{ii}]; end
    url = ['https://raw.githubusercontent.com/uwa-channels/matlab/' revision '/' relative];
    if ~isfile(filename)
        websave([filename '.part'], url, weboptions('Timeout',60));
        movefile([filename '.part'], filename);
    end
    if ~isempty(previous_code)
        old = find(strcmp({previous_code.name},code_names{ii}));
        assert(isscalar(old) && strcmp(black_file_hash(filename),previous_code(old).sha256), ...
            'black:SourceDigest','Pinned local source changed: %s',code_names{ii});
    end
    manifest.code(ii) = struct('name',code_names{ii},'url',url,'sha256',black_file_hash(filename));
end
fid = fopen(manifest_file,'w','n','UTF-8');
assert(fid >= 0, 'black:ManifestWrite', 'Cannot write manifest.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(manifest,PrettyPrint=true));
assets = struct('case_dir',case_dir,'raw_dir',raw_dir,'vendor_dir',vendor_dir, ...
    'manifest_file',manifest_file,'manifest',manifest);
end
