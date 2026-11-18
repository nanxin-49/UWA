function receipt=audit_r2_seed260002_saved(source,out,mode)
%AUDIT_R2_SEED260002_SAVED Read-only retrieval diagnosis, never solve/restart.
arguments
    source (1,:) char
    out (1,:) char
    mode (1,:) char = 'numeric'
end
assert(ismember(mode,{'numeric','figures','renderer-replay'}));
assert(isfile(fullfile(out,'integrity_receipt.json')));
integrity=jsondecode(fileread(fullfile(out,'integrity_receipt.json')));
assert(strcmp(integrity.status,'LOCAL_RETRIEVAL_AND_CHAIN_SHA_PASS'));
source=char(java.io.File(source).getCanonicalPath()); out=char(java.io.File(out).getCanonicalPath());
assert(~startsWith(out,[source filesep]));
oldpath=path; pathguard=onCleanup(@()path(oldpath));
threads=maxNumCompThreads(1); threadguard=onCleanup(@()maxNumCompThreads(threads));
addpath(fullfile(source,'code'),'-begin');
guards=fullfile(out,'solve_sentinels'); if ~isfolder(guards), mkdir(guards); end
for name={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
        'closeout_s3a_pm_fullfk','run_pe_1d_surface_reflection_validation','run_r2_piston_pm_job'}
    p=fullfile(guards,[name{1} '.m']);
    if ~isfile(p)
        fid=fopen(p,'w'); assert(fid>=0);
        fprintf(fid,'function varargout=%s(varargin)\nerror(''r2review:ForbiddenSolve'',''Saved-data review only'');\nend\n',name{1}); fclose(fid);
    end
end
addpath(guards,'-begin'); rehash;
assert(startsWith(which('solve_helmholtz_bie_halfplane_vertical'),guards));
timer=tic; diary(fullfile(out,[mode '_diary.log'])); dg=onCleanup(@()diary('off'));
work=fullfile(source,'finalization','attempts','attempt_0001','work');
if strcmp(mode,'numeric')
    dest=fullfile(out,'reaudit'); assert(~isfolder(dest)); mkdir(dest);
    inventory=whos('-file',fullfile(work,'common.mat'));
    assert(all(ismember({inventory.class},{'struct','table','double','cell','logical'})));
    for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','policy.json'}
        copyfile(fullfile(work,name{1}),fullfile(dest,name{1}));
    end
    ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
    for j=1:5
        d=fullfile(dest,'jobs',ids{j}); mkdir(d);
        copyfile(fullfile(source,'jobs',ids{j},'data.mat'),fullfile(d,'data.mat'));
        copyfile(fullfile(source,'jobs',ids{j},'points.csv'),fullfile(d,'points.csv'));
    end
    c3=audit_r1_bie_manual_stage(fullfile(dest,'jobs','bie_space'),fullfile(dest,'common.mat'),[60 14 96]);
    c4=audit_r1_bie_manual_stage(fullfile(dest,'jobs','bie_window'),fullfile(dest,'common.mat'),[70 12 96]);
    comprehensive=audit_piston_pm_r1(dest); diary(fullfile(out,'numeric_diary.log'));
    final=audit_r2_pm_final(dest);
    summary=summarize_r2_pm_fields({dest},fullfile(out,'seed260002_spatial_summary'));
    assert(c3.failed_checks==0 && c4.failed_checks==0 && comprehensive.failed_checks==0 && final.failed_checks==0);
    receipt=struct('status','SAVED_NUMERICAL_REAUDIT_PASS','failed_checks',0,'elapsed_s',toc(timer), ...
        'MATLAB_version',version,'C3',c3,'C4',c4,'comprehensive',comprehensive,'final',final, ...
        'summary_status',summary.status,'new_density_PS_FullFK_solves',0,'SSH_connections',0, ...
        'scope','C3/C4 independent saved density audits, comprehensive source/kernel/CSV, four-pair resolution and spatial metrics; no new solve');
elseif strcmp(mode,'figures')
    dest=fullfile(out,'figure_diagnostics'); assert(~isfolder(dest)); mkdir(dest);
    q=load(fullfile(work,'common.mat'),'x'); a=load(fullfile(work,'comparison.mat'));
    sourcefig=fullfile(work,'figures_final'); records=struct; total=0;
    for item={'F02_pointwise_errors','F04_FullFK_pointwise_errors'}
        name=item{1}; f=openfig(fullfile(sourcefig,[name '.fig']),'new','invisible'); cg=onCleanup(@()close(f));
        lines=findall(f,'Type','line'); count=0;
        for j=1:numel(lines)
            info=lines(j).UserData; if ~isstruct(info) || ~isfield(info,'panel'), continue; end
            m=3; if strcmp(info.method,'PS'), m=1; end
            e=a.U(:,m)-a.U(:,2);
            switch info.panel
                case 1, v=abs(e);
                case 2, v=abs(e)/a.scale;
                case 3, v=abs(a.U(:,m))-abs(a.U(:,2));
                case 4, v=angle(a.U(:,m).*conj(a.U(:,2))); v(~a.mask)=NaN;
            end
            assert(isequaln(lines(j).XData(:),q.x) && isequaln(lines(j).YData(:),v)); count=count+1;
        end
        expected=4; if startsWith(name,'F02'), expected=8; end
        assert(count==expected); total=total+count;
        before=local_state(f); f.Position=[80 80 1150 1000]; drawnow; after_draw=local_state(f);
        states=cell(3,1); exports=cell(3,1);
        for j=1:3
            drawnow; stem=sprintf('%s_probe%d',name,j);
            savefig(f,fullfile(dest,[stem '.fig']));
            exportgraphics(f,fullfile(dest,[stem '.png']),'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
            states{j}=local_state(f); exports{j}=imread(fullfile(dest,[stem '.png']));
        end
        clear cg;
        original=imread(fullfile(sourcefig,[name '.png']));
        r=struct('checked_error_curves',count,'before_draw',before,'after_draw',after_draw,'after_exports',{states}, ...
            'original_vs_local_export1',local_pixels(original,exports{1}), ...
            'local_export1_vs_export2',local_pixels(exports{1},exports{2}), ...
            'local_export2_vs_export3',local_pixels(exports{2},exports{3}));
        if startsWith(name,'F02')
            old=imread(fullfile(sourcefig,'F02_audit_reexport.png'));
            r.original_remote_pair=local_pixels(original,old);
            r.remote_reexport_vs_local_export1=local_pixels(old,exports{1});
        end
        records.(name)=r;
    end
    receipt=struct('status','FIG_DATA_AND_EXPORT_DIAGNOSIS_COMPLETE','failed_data_checks',0,'checked_error_curves',total, ...
        'elapsed_s',toc(timer),'records',records,'new_acoustic_solves',0,'SSH_connections',0, ...
        'original_failure_preserved',true,'full_seed_final_seal_created',false);
else
    dest=fullfile(out,'reaudit'); name='original_renderer_replay';
    assert(~isfolder(fullfile(dest,name)));
    try
        r=render_r2_pm_seed(dest,name);
        receipt=struct('status','ORIGINAL_RENDERER_REPLAY_PASS','renderer_receipt',r,'elapsed_s',toc(timer));
    catch ME
        receipt=struct('status','ORIGINAL_RENDERER_FAILURE_REPRODUCED','identifier',ME.identifier, ...
            'message',ME.message,'stack',ME.stack,'elapsed_s',toc(timer));
        fprintf('%s\n',getReport(ME,'extended','hyperlinks','off'));
    end
    receipt.new_acoustic_solves=0; receipt.SSH_connections=0;
end
r1_bie_manual_io('write',fullfile(out,[mode '_receipt.json']),receipt);
fprintf('REVIEW %s %.3fs\n',receipt.status,toc(timer));
end
function s=local_state(f)
ax=findall(f,'Type','axes'); positions=zeros(numel(ax),4); xl=zeros(numel(ax),2); yl=xl;
for j=1:numel(ax), positions(j,:)=ax(j).Position; xl(j,:)=ax(j).XLim; yl(j,:)=ax(j).YLim; end
s=struct('figure_units',f.Units,'figure_position',f.Position,'axes_positions',positions,'xlimits',xl,'ylimits',yl);
end
function s=local_pixels(a,b)
s=struct('same_dimensions',isequal(size(a),size(b)),'dimensions_a',size(a),'dimensions_b',size(b));
if ~s.same_dimensions, return; end
d=any(a~=b,3); [y,x]=find(d); bbox=[]; if ~isempty(x), bbox=[min(x) min(y) max(x) max(y)]; end
s.changed_pixels=nnz(d); s.changed_channel_values=nnz(a~=b); s.changed_fraction=nnz(d)/numel(d);
s.maximum_channel_difference=max(abs(double(a)-double(b)),[],'all'); s.changed_bbox_xy=bbox;
end
