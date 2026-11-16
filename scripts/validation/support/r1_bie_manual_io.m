function out=r1_bie_manual_io(action,varargin)
%R1_BIE_MANUAL_IO Hash-bound, atomic checkpoint primitives; no physics.
switch action
    case 'sha'
        path=varargin{1}; fid=fopen(path,'rb'); assert(fid>=0,'Cannot read %s',path);
        guard=onCleanup(@()fclose(fid));
        md=java.security.MessageDigest.getInstance('SHA-256');
        while ~feof(fid), bytes=fread(fid,1048576,'*uint8'); md.update(typecast(bytes,'int8')); end
        out=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
    case 'read'
        out=jsondecode(fileread(varargin{1}));
    case 'write'
        path=varargin{1}; value=varargin{2}; temp=[path '.tmp'];
        fid=fopen(temp,'w','n','UTF-8'); assert(fid>=0,'Cannot write %s',temp);
        fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true)); fclose(fid);
        [ok,msg]=movefile(temp,path,'f'); assert(ok,'%s',msg); out=path;
    case 'path'
        root=varargin{1}; relative=strrep(char(varargin{2}),'\','/');
        assert(~isempty(relative) && isempty(regexp(relative,'(^/|:|(^|/)\.\.(/|$))','once')), ...
            'r1manual:UnsafePath','Invalid relative artifact path');
        out=fullfile(root,strrep(relative,'/',filesep));
    case 'verify_manifest'
        root=varargin{1}; out=r1_bie_manual_io('read',fullfile(root,'manifest.json'));
        assert(out.schema_version==1 && strcmp(out.scope,'R1_BIE_C1_C4_MANUAL'));
        for j=1:numel(out.files)
            f=out.files(j); path=r1_bie_manual_io('path',root,f.path);
            assert(strcmp(r1_bie_manual_io('sha',path),f.sha256), ...
                'r1manual:InputDrift','Package input/code SHA mismatch: %s',f.path);
        end
    case 'seal'
        root=varargin{1}; name=varargin{2}; signature=varargin{3}; paths=varargin{4}; status=varargin{5};
        target=fullfile(root,name); assert(~isfile(target),'Refusing to replace sealed %s',target);
        files=struct('path',{},'sha256',{});
        for j=1:numel(paths)
            rel=strrep(paths{j},'\','/'); path=r1_bie_manual_io('path',root,rel);
            files(end+1)=struct('path',rel,'sha256',r1_bie_manual_io('sha',path)); %#ok<AGROW>
        end
        out=struct('schema_version',1,'status',status,'signature',signature,'files',files);
        r1_bie_manual_io('write',target,out);
    case 'verify_seal'
        root=varargin{1}; name=varargin{2}; signature=varargin{3};
        out=r1_bie_manual_io('read',fullfile(root,name));
        assert(out.schema_version==1 && strcmp(out.signature,signature), ...
            'r1manual:StaleCheckpoint','Checkpoint inputs/code/predecessor changed');
        assert(~isempty(out.files),'Empty checkpoint artifact binding');
        for j=1:numel(out.files)
            f=out.files(j); path=r1_bie_manual_io('path',root,f.path);
            assert(strcmp(r1_bie_manual_io('sha',path),f.sha256), ...
                'r1manual:ArtifactDrift','Checkpoint artifact mismatch: %s',f.path);
        end
    case 'resume_state'
        root=varargin{1}; signature=varargin{2};
        if isfile(fullfile(root,'committed.json'))
            seal=r1_bie_manual_io('verify_seal',root,'committed.json',signature);
            assert(strcmp(seal.status,'AUDITED_PASS'));
            audit=r1_bie_manual_io('read',fullfile(root,'audit.json'));
            gates=r1_bie_manual_io('read',fullfile(root,'gates.json'));
            assert(strcmp(audit.status,'PASS') && audit.failed_checks==0 && strcmp(gates.status,'PASS'));
            out='SKIP_AUDITED';
        elseif isfile(fullfile(root,'solved.json'))
            seal=r1_bie_manual_io('verify_seal',root,'solved.json',signature);
            assert(strcmp(seal.status,'SOLVED_UNAUDITED')); out='AUDIT_ONLY';
        elseif isfile(fullfile(root,'data.mat')) || isfile(fullfile(root,'points.csv'))
            error('r1manual:OrphanData','Unsealed data exists; stop for review, never silently recompute: %s',root);
        else
            out='SOLVE';
        end
    otherwise
        error('Unknown manual checkpoint action: %s',action);
end
end
