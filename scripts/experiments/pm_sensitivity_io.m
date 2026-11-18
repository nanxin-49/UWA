function out = pm_sensitivity_io(action, varargin)
%PM_SENSITIVITY_IO Research-only portable provenance and checkpoint operations.
switch action
    case 'sha'
        fid=fopen(varargin{1},'rb'); assert(fid>=0); guard=onCleanup(@()fclose(fid));
        md=java.security.MessageDigest.getInstance('SHA-256');
        while ~feof(fid)
            bytes=fread(fid,1048576,'*uint8');
            if ~isempty(bytes), md.update(typecast(bytes,'int8')); end
        end
        out=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
    case 'write'
        file=varargin{1}; tmp=[file '.pending'];
        fid=fopen(tmp,'w','n','UTF-8'); assert(fid>=0);
        guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(varargin{2},PrettyPrint=true));
        clear guard; movefile(tmp,file,'f'); out=file;
    case 'log'
        fid=fopen(varargin{1},'a','n','UTF-8'); assert(fid>=0);
        guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',varargin{2});
        fprintf('%s\n',varargin{2}); out=varargin{1};
    case 'manifest'
        root=varargin{1}; paths=varargin{2};
        out=struct('path',{},'sha256',{});
        for j=1:numel(paths)
            out(j)=struct('path',strrep(paths{j},'\','/'), ...
                'sha256',pm_sensitivity_io('sha',fullfile(root,paths{j})));
        end
    case 'verify'
        root=varargin{1}; list=varargin{2};
        for j=1:numel(list)
            p=list(j).path;
            assert(isempty(regexp(p,'(^/|:|\\|(^|/)\.\.(/|$))','once')),'Unsafe relative path');
            assert(strcmp(pm_sensitivity_io('sha',fullfile(root,p)),list(j).sha256),'SHA drift: %s',p);
        end
        out=true;
    case 'absolute'
        out=varargin{1}; assert(~isempty(regexp(out,'^([A-Za-z]:[\\/]|/)','once')),'Absolute path required');
    otherwise
        error('Unknown IO operation');
end
end
