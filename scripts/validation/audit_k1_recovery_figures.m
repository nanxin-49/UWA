function receipt=audit_k1_recovery_figures(work_dirs,figdir,out)
%AUDIT_K1_RECOVERY_FIGURES Independently verify 12 saved curves and zero pixels.
arguments
    work_dirs (1,1) cell
    figdir (1,:) char
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out);
q=load(fullfile(work_dirs{1},'common.mat'),'cfg','x'); a=load(fullfile(work_dirs{1},'comparison.mat'));
assert(a.reference_qualified && numel(q.x)==411);
assert(abs(a.scale-norm(a.U(:,2))/sqrt(411))<1e-12 && isequal(a.mask,abs(a.U(:,2))>=.01*a.scale));
kinds={'absolute','normalized','amplitude','phase'}; records=struct; total=0;
try
for family=1:2
    prefix='F02'; expected=2; if family==2, prefix='F04'; expected=1; end
    for panel=1:4
        stem=[prefix '_' kinds{panel}]; original=imread(fullfile(figdir,[stem '.png']));
        assert(isequal(size(original),[550 1000 3])); layouts=cell(1,2); differences=zeros(1,2);
        for pass=1:2
            f=openfig(fullfile(figdir,[stem '.fig']),'new','invisible'); guard=onCleanup(@()local_close(f));
            f.Units='pixels'; f.Position=[10 10 1000 550]; drawnow;
            lines=findall(f,'Type','line'); assert(numel(lines)==expected); methods={};
            for j=1:numel(lines)
                u=lines(j).UserData; assert(u.seed==q.cfg.seed && u.panel==panel && strcmp(u.kind,kinds{panel}));
                m=1; if strcmp(u.method,'FullFK'), m=3; else, assert(strcmp(u.method,'PS')); end
                if family==2, assert(m==3); end
                b=a.U(:,2); v=a.U(:,m);
                switch kinds{panel}
                    case 'absolute', y=abs(v-b);
                    case 'normalized', y=abs(v-b)/a.scale;
                    case 'amplitude', y=abs(v)-abs(b);
                    case 'phase', y=angle(v.*conj(b)); y(abs(b)<.01*a.scale)=NaN;
                end
                assert(isequaln(lines(j).XData(:),q.x) && isequaln(lines(j).YData(:),y));
                methods{end+1}=u.method; %#ok<AGROW>
            end
            assert(numel(unique(methods))==expected);
            before=local_layout(f,q.x); target=fullfile(out,sprintf('%s_reopen%d.png',stem,pass));
            f.PaperUnits='inches'; f.PaperPosition=[0 0 10 5.5]; f.PaperSize=[10 5.5]; f.PaperPositionMode='manual';
            print(f,target,'-dpng','-r100');
            assert(isequaln(before,local_layout(f,q.x)),'k1small:Layout','Layout changed during export');
            layouts{pass}=before; rendered=imread(target); assert(isequal(size(original),size(rendered)));
            differences(pass)=nnz(original~=rendered); records.(stem)=struct('pixel_differences',differences,'curves',expected);
            assert(differences(pass)==0,'k1small:Pixels','Nonzero PNG channel differences: %d',differences(pass)); clear guard;
        end
        assert(isequaln(layouts{1},layouts{2})); total=total+expected; records.(stem).layout_exact=true;
    end
end
assert(total==12);
receipt=struct('status','K1_SMALL_PANEL_FIG_AND_ZERO_PIXEL_PASS','failed_checks',0,'layout_version','K1-small-panels-v3-fixed-paper','checked_error_curves',total,'records',records,'new_acoustic_solves',0,'pixel_tolerance',0);
r1_bie_manual_io('write',fullfile(out,'receipt.json'),receipt);
catch ME
    r1_bie_manual_io('write',fullfile(out,'failure.json'),struct('status','FAIL','error',ME.message,'records',records)); rethrow(ME);
end
end
function s=local_layout(f,x)
ax=findall(f,'Type','axes'); lg=findall(f,'Type','legend'); assert(isscalar(ax) && isscalar(lg));
assert(isequal(f.Position,[10 10 1000 550]) && isequal(ax.Position,[100 75 700 400]) && isequal(lg.Position,[805 365 140 75]) && ax.FontSize==11 && lg.FontSize==10);
assert(isequal(ax.XLim,[min(x) max(x)]) && isequal(ax.XTick,-40:10:40));
assert(strcmp(ax.XLimMode,'manual') && strcmp(ax.YLimMode,'manual') && strcmp(ax.XTickMode,'manual') && strcmp(ax.YTickMode,'manual') && strcmp(ax.XTickLabelMode,'manual') && strcmp(ax.YTickLabelMode,'manual'));
s=struct('figure',f.Position,'axes',ax.Position,'legend',lg.Position,'limits',[ax.XLim ax.YLim],'fonts',[ax.FontSize lg.FontSize]);
end
function local_close(f)
if isgraphics(f), close(f); end
end
