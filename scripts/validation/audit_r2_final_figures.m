function receipt=audit_r2_final_figures(work_dirs,figdir,out,mode)
%AUDIT_R2_FINAL_FIGURES Independent saved-array checks and strict zero pixels.
arguments
    work_dirs (1,:) cell
    figdir (1,:) char
    out (1,:) char
    mode (1,:) char = 'stable'
end
assert(~isfolder(out) && ismember(mode,{'stable','original-data'})); mkdir(out);
Q=cell(size(work_dirs)); A=Q;
for j=1:numel(Q), Q{j}=load(fullfile(work_dirs{j},'common.mat'),'cfg','x'); A{j}=load(fullfile(work_dirs{j},'comparison.mat')); end
names={'F02_pointwise_errors','F04_FullFK_pointwise_errors'}; expected=[8 4];
if numel(Q)==3, names={'F05_three_seed_errors'}; expected=12; end
records=struct; total=0;
try
for k=1:numel(names)
    stem=names{k}; f=openfig(fullfile(figdir,[stem '.fig']),'new','invisible'); guard=onCleanup(@()local_close(f));
    lines=findall(f,'Type','line'); keys={};
    for j=1:numel(lines)
        u=lines(j).UserData; if ~isstruct(u) || ~isfield(u,'panel'), continue; end
        if ~isfield(u,'seed'), u.seed=Q{1}.cfg.seed; end
        ix=find(cellfun(@(q)q.cfg.seed==u.seed,Q)); assert(isscalar(ix));
        if ~isfield(u,'kind'), kinds={'absolute','normalized','amplitude','phase'}; u.kind=kinds{u.panel}; end
        m=3; if strcmp(u.method,'PS'), m=1; else, assert(strcmp(u.method,'FullFK')); end
        b=A{ix}.U(:,2); model=A{ix}.U(:,m);
        scale=A{ix}.scale;
        assert(abs(scale-norm(b)/sqrt(411))<1e-12 && isequal(A{ix}.mask,abs(b)>=.01*scale));
        switch u.kind
            case 'absolute', v=abs(model-b);
            case 'normalized', v=abs(model-b)/scale;
            case 'amplitude', v=abs(model)-abs(b);
            case 'phase', v=angle(model.*conj(b)); v(~A{ix}.mask)=NaN;
            otherwise, error('r2figure:Kind','Unknown curve semantic');
        end
        assert(isequaln(lines(j).XData(:),Q{ix}.x) && isequaln(lines(j).YData(:),v),'r2figure:Data','FIG curve mismatch');
        keys{end+1}=sprintf('%d:%s:%d',u.seed,u.method,u.panel); %#ok<AGROW>
    end
    assert(numel(keys)==expected(k) && numel(unique(keys))==expected(k)); total=total+numel(keys);
    rec=struct('checked_curves',numel(keys),'keys',{keys});
    if strcmp(mode,'stable')
        rec.position_after_openfig=f.Position;
        % openfig may shrink even invisible figures to fit the host screen.
        % Restore the declared raster canvas, not any scientific array.
        f.Units='pixels'; f.Position=[40 40 1800 1400]; drawnow;
        before=local_layout(f); original=imread(fullfile(figdir,[stem '.png'])); pixel_values=zeros(1,2);
        for pass=1:2
            if pass==2
                clear guard; f=openfig(fullfile(figdir,[stem '.fig']),'new','invisible'); guard=onCleanup(@()local_close(f));
                f.Units='pixels'; f.Position=[40 40 1800 1400];
            end
            drawnow; audit_stem=sprintf('%s_reopen%d',stem,pass);
            savefig(f,fullfile(out,[audit_stem '.fig']));
            exportgraphics(f,fullfile(out,[audit_stem '.png']),'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
            rendered=imread(fullfile(out,[audit_stem '.png'])); assert(isequal(size(original),size(rendered)));
            pixel_values(pass)=nnz(original~=rendered);
            rec.pixel_differences=pixel_values; records.(stem)=rec;
            assert(isequaln(before,local_layout(f)),'r2figure:Layout','Figure/axes layout changed during export');
            assert(pixel_values(pass)==0,'r2figure:Pixels','Nonzero PNG pixel differences: %d',pixel_values(pass));
        end
        rec.pixel_differences=pixel_values; rec.layout_exact=true;
    elseif k==1 && isfile(fullfile(figdir,'F02_audit_reexport.png'))
        im1=imread(fullfile(figdir,[stem '.png'])); im2=imread(fullfile(figdir,'F02_audit_reexport.png'));
        assert(isequal(size(im1),size(im2))); rec.historical_pixel_differences=nnz(im1~=im2);
        rec.historical_failure_not_overridden=true;
    end
    records.(stem)=rec; clear guard;
end
receipt=struct('status','FIG_ARRAYS_AND_ZERO_PIXEL_AUDIT_PASS','failed_checks',0,'mode',mode, ...
    'checked_error_curves',total,'records',records,'new_acoustic_solves',0,'pixel_tolerance',0);
if strcmp(mode,'original-data'), receipt.status='ORIGINAL_FIG_DATA_PASS_PIXEL_HISTORY_PRESERVED'; end
r1_bie_manual_io('write',fullfile(out,'receipt.json'),receipt);
catch ME
    r1_bie_manual_io('write',fullfile(out,'failure.json'),struct('status','FAIL','error',ME.message,'identifier',ME.identifier,'records',records));
    rethrow(ME);
end
end
function s=local_layout(f)
ax=findall(f,'Type','axes'); positions=zeros(numel(ax),4); limits=zeros(numel(ax),4);
for j=1:numel(ax)
    assert(strcmp(ax(j).Units,'pixels') && strcmp(ax(j).XLimMode,'manual') && strcmp(ax(j).YLimMode,'manual') && ...
        strcmp(ax(j).XTickMode,'manual') && strcmp(ax(j).YTickMode,'manual'));
    positions(j,:)=ax(j).Position; limits(j,:)=[ax(j).XLim ax(j).YLim];
end
s=struct('figure_position',f.Position,'positions',positions,'limits',limits);
end
function local_close(f)
if isgraphics(f), close(f); end
end
