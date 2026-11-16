function out=r1_bie_continuation_io(action,varargin)
%R1_BIE_CONTINUATION_IO Separate v2 provenance chain; never rewrites v1.
switch action
    case 'sha'
        out=local_sha(varargin{1});
    case {'read','write','path'}
        out=r1_bie_manual_io(action,varargin{:});
    case 'policy'
        out=varargin{1};
        assert(out.schema_version==2 && strcmp(out.version,'r1-bie-reference-v2') && ...
            strcmp(out.grid_monotonicity_mode,'diagnostic') && out.reference_stage==3);
        expected=struct('sensitivity',1e-4,'linear',1e-10,'boundary',1e-4, ...
            'reconstruction',1e-11,'geometry',1e-12,'source',1e-7, ...
            'PS',1e-4,'PS_step',1e-10,'FullFK',1e-5,'resolution_factor',5,'phase_mask_fraction',.01);
        assert(isequaln(orderfields(out.gates),orderfields(expected)), ...
            'r1v2:Policy','The approved numerical gates must not change');
        assert(isequal(string(out.pair_names(:)),["grid_10_12";"grid_12_14";"self_64_96";"window_60_70"]));
    case 'verify_package'
        root=varargin{1}; out=r1_bie_manual_io('read',fullfile(root,'manifest.json'));
        assert(out.schema_version==2 && strcmp(out.scope,'R1_BIE_C4_CONTINUATION'));
        assert(out.allowed_new_solve_stage==4 && out.single_process && out.single_thread && ...
            ~out.parpool && out.RAM_reserve_fraction==.3 && ~out.automatic_retry);
        local_files(root,out.files);
        policy=r1_bie_continuation_io('policy',r1_bie_manual_io('read',fullfile(root,'policy.json')));
        assert(strcmp(out.policy_sha256,r1_bie_continuation_io('sha',fullfile(root,'policy.json'))));
        approval=r1_bie_manual_io('read',fullfile(root,'approval.json'));
        assert(approval.implementation_approved && ~approval.C4_started_by_agent && ...
            strcmp(approval.policy_version,policy.version));
    case 'verify_origin'
        root=varargin{1}; origin=fullfile(root,'origin_evidence');
        sources=r1_bie_manual_io('read',fullfile(root,'migration_sources.json'));
        assert(strcmp(r1_bie_continuation_io('sha',fullfile(origin,'manifest.json')),sources.legacy_manifest_sha256));
        list=r1_bie_manual_io('path',origin,sources.retrieval_manifest_path);
        assert(strcmp(r1_bie_continuation_io('sha',list),sources.retrieval_manifest_sha256));
        lines=splitlines(strtrim(fileread(list))); paths=cell(numel(lines),1);
        for j=1:numel(lines)
            row=regexp(char(lines(j)),'^([0-9a-f]{64})  (.+)$','tokens','once');
            assert(~isempty(row),'Malformed legacy retrieval SHA'); paths{j}=row{2};
            assert(strcmp(r1_bie_continuation_io('sha',r1_bie_manual_io('path',origin,row{2})),row{1}), ...
                'r1v2:OriginDrift','Original retrieval artifact drift: %s',row{2});
        end
        assert(numel(paths)==110 && numel(unique(paths))==110);
        legacy=r1_bie_manual_io('read',fullfile(origin,'manifest.json'));
        assert(legacy.schema_version==1 && strcmp(legacy.scope,'R1_BIE_C1_C4_MANUAL') && numel(legacy.files)==36); local_files(origin,legacy.files);
        pred=sources.legacy_manifest_sha256; ids=local_ids();
        for j=1:4
            dest=fullfile(origin,'jobs',ids{j}); sig=[pred ':C' num2str(j-1) ':' sources.legacy_manifest_sha256];
            solved=local_legacy_seal(dest,'solved.json',sig);
            assert(strcmp(solved.status,'SOLVED_UNAUDITED') && numel(solved.files)==4);
            local_required(solved.files,{'data.mat','points.csv','config.json','timings.json'});
            audit=r1_bie_manual_io('read',fullfile(dest,'audit.json'));
            assert(strcmp(audit.status,'PASS') && audit.failed_checks==0 && ...
                numel(fieldnames(audit.checks))==6 && all(cell2mat(struct2cell(audit.checks))));
            required={'data.mat','points.csv','config.json','timings.json','audit.json','audit_reconstruction.mat','gates.json','solved.json'};
            for k=1:numel(required)
                target=['jobs/' ids{j} '/' required{k}];
                assert(any(strcmp(paths,target)) || any(strcmp(paths,['./' target])), ...
                    'r1v2:UnboundEvidence','Legacy evidence not in retrieval SHA: %s',target);
            end
            gates=r1_bie_manual_io('read',fullfile(dest,'gates.json'));
            if j<4
                seal=local_legacy_seal(dest,'committed.json',sig);
                local_required(seal.files,required(1:7));
                assert(strcmp(seal.status,'AUDITED_PASS') && strcmp(gates.status,'PASS'));
                pred=r1_bie_continuation_io('sha',fullfile(dest,'committed.json'));
            else
                assert(~isfile(fullfile(dest,'committed.json')) && strcmp(gates.status,'FAIL') && gates.grid_nonmonotonic);
                assert(numel(gates.pairs)==3 && all([gates.pairs.relative]<1e-4));
            end
        end
        cp=r1_bie_manual_io('read',fullfile(origin,'checkpoint.json'));
        assert(strcmp(cp.error_identifier,'r1manual:Nonmonotonic') && cp.new_acoustic_solves==3);
        assert(~isfile(fullfile(origin,'jobs','bie_window','data.mat')));
        out=struct('status','PASS','retrieval_SHA_files',110,'legacy_package_SHA_files',36, ...
            'legacy_policy_status','FAIL','historical_completed_solves',3);
    case 'seal'
        root=varargin{1}; name=varargin{2}; sig=varargin{3}; paths=varargin{4}; status=varargin{5};
        target=fullfile(root,name); assert(~isfile(target),'r1v2:ReplaceSeal','Refusing to replace a v2 seal');
        files=struct('path',{},'sha256',{});
        for j=1:numel(paths)
            rel=strrep(paths{j},'\','/');
            files(end+1)=struct('path',rel,'sha256',r1_bie_continuation_io('sha',r1_bie_manual_io('path',root,rel))); %#ok<AGROW>
        end
        out=struct('schema_version',2,'status',status,'signature',sig,'files',files);
        r1_bie_manual_io('write',target,out);
    case 'verify_seal'
        root=varargin{1}; name=varargin{2}; sig=varargin{3}; required=varargin{4};
        out=r1_bie_manual_io('read',fullfile(root,name));
        assert(out.schema_version==2 && strcmp(out.signature,sig),'r1v2:Stale','v2 inputs/policy/predecessor changed');
        local_required(out.files,required); local_files(root,out.files);
    case 'verify_migration'
        root=varargin{1}; r1_bie_continuation_io('verify_package',root);
        r1_bie_continuation_io('verify_origin',root);
        fp=r1_bie_continuation_io('sha',fullfile(root,'manifest.json'));
        ph=r1_bie_continuation_io('sha',fullfile(root,'policy.json')); pred=fp; ids=local_ids();
        required={'data.mat','points.csv','config.json','timings.json','audit.json', ...
            'audit_reconstruction.mat','gates_v2.json','effective_config.json'};
        for j=1:4
            dest=fullfile(root,'jobs',ids{j}); sig=[pred ':C' num2str(j-1) ':' fp ':' ph];
            seal=r1_bie_continuation_io('verify_seal',dest,'migration_commit.json',sig,required);
            expected='IMPORTED_AUDITED_PASS';
            if j==4, expected='IMPORTED_AUDITED_PASS_WITH_NONMONOTONIC_DIAGNOSTIC'; end
            assert(strcmp(seal.status,expected));
            for name={'data.mat','points.csv','config.json','timings.json'}
                assert(strcmp(r1_bie_continuation_io('sha',fullfile(dest,name{1})), ...
                    r1_bie_continuation_io('sha',fullfile(root,'origin_evidence','jobs',ids{j},name{1}))));
            end
            a=r1_bie_manual_io('read',fullfile(dest,'audit.json'));
            g=r1_bie_manual_io('read',fullfile(dest,'gates_v2.json'));
            assert(strcmp(a.status,'PASS') && a.failed_checks==0 && strcmp(g.status,'PASS') && strcmp(g.policy_sha256,ph));
            pred=r1_bie_continuation_io('sha',fullfile(dest,'migration_commit.json'));
        end
        marker=r1_bie_continuation_io('verify_seal',root,'migration_complete.json',[pred ':MIGRATED:' fp ':' ph], ...
            {'migration_checkpoint.json','migration_audits/C3/audit.json','migration_audits/C3/audit_reconstruction.mat'});
        assert(strcmp(marker.status,'MIGRATION_AUDITED_PASS'));
        out=struct('predecessor',pred,'manifest_sha256',fp,'policy_sha256',ph);
    case 'c4_state'
        root=varargin{1}; sig=varargin{2};
        payload={'data.mat','points.csv','config.json','timings.json','solve_invocation.json'};
        if isfile(fullfile(root,'committed.json'))
            seal=r1_bie_continuation_io('verify_seal',root,'committed.json',sig, ...
                [payload {'solved.json','audit.json','audit_reconstruction.mat','gates_v2.json'}]);
            a=r1_bie_manual_io('read',fullfile(root,'audit.json')); g=r1_bie_manual_io('read',fullfile(root,'gates_v2.json'));
            assert(strcmp(seal.status,'AUDITED_PASS') && strcmp(a.status,'PASS') && a.failed_checks==0 && strcmp(g.status,'PASS'));
            out='SKIP_AUDITED';
        elseif isfile(fullfile(root,'solved.json'))
            seal=r1_bie_continuation_io('verify_seal',root,'solved.json',sig,payload);
            assert(strcmp(seal.status,'SOLVED_UNAUDITED')); out='AUDIT_ONLY';
        elseif isfile(fullfile(root,'solve_invocation.json')) || isfile(fullfile(root,'data.mat')) || isfile(fullfile(root,'points.csv'))
            error('r1v2:OrphanData','C4 invoked or unsealed data present: stop, never automatically solve again');
        else
            out='SOLVE_C4_ONCE';
        end
    otherwise
        error('Unknown v2 IO action: %s',action);
end
end
function local_files(root,files)
assert(~isempty(files) && numel(unique(string({files.path})))==numel(files));
for j=1:numel(files)
    f=files(j); assert(strcmp(r1_bie_continuation_io('sha',r1_bie_manual_io('path',root,f.path)),f.sha256), ...
        'r1v2:ArtifactDrift','SHA mismatch: %s',f.path);
end
end
function local_required(files,names)
assert(all(ismember(string(names),string({files.path}))),'r1v2:MissingBinding','Required evidence missing from seal');
end
function ids=local_ids()
ids={'bie_coarse','bie_main','bie_self','bie_space'};
end
function out=local_sha(path)
fid=fopen(path,'rb'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while true
    bytes=fread(fid,1048576,'*uint8');
    if isempty(bytes), break; end
    md.update(typecast(bytes,'int8'));
end
out=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function out=local_legacy_seal(root,name,signature)
out=r1_bie_manual_io('read',fullfile(root,name));
assert(out.schema_version==1 && strcmp(out.signature,signature));
local_files(root,out.files);
end
