function summary = generate_bellhop_fig12b_three_attempt_visuals(outdir)
%GENERATE_BELLHOP_FIG12B_THREE_ATTEMPT_VISUALS Raw saved-field comparison.
arguments
    outdir (1,:) char = ''
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(outdir), outdir=fullfile(root,'results','validation','bellhop_fig12b_three_attempts'); end
figdir=fullfile(outdir,'figures'); if ~isfolder(figdir), mkdir(figdir); end
labels={'Geometric Gaussian / point source','Cerveny ray-centered / point source', ...
    'Geometric Gaussian / line source'};
summary=struct([]); checks=struct([]); controls=struct([]);
fig=figure('Visible','off','Color','w','Position',[40 40 1700 1680]);
layout=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
for attempt=1:3
    tags=arrayfun(@(n)sprintf('a%d_n%d_h1.25',attempt,n),[50 500 2000],'UniformOutput',false);
    tags{4}=sprintf('a%d_n2000_h0.625',attempt); fields=cell(1,4);
    for jj=1:4
        s=load(fullfile(outdir,[tags{jj} '.mat']),'output'); o=s.output; fields{jj}=o;
        ax=nexttile(layout); local_field(ax,o);
        title(ax,sprintf('A%d: N=%d, step=%g m',attempt,o.cfg.beam_count,o.cfg.step_m),'Interpreter','none');
        mask=o.depth_m>=10 & o.depth_m<=4900 & abs(o.range_m)>=.1;
        background=o.depth_m>=10 & o.depth_m<=3500 & abs(o.range_m)>=.5;
        tl=-20*log10(abs(o.pressure)); values=tl(background & isfinite(tl));
        summary=[summary struct('attempt',attempt,'tag',string(tags{jj}), ...
            'run_type',string(o.cfg.run_type),'beams',o.cfg.beam_count,'step_m',o.cfg.step_m, ...
            'seconds',o.elapsed_s,'zero_fraction',o.zero_fraction,'nonfinite_count',o.nonfinite_count, ...
            'storage_warnings',o.storage_warnings,'unbounded_warnings',o.unbounded_warnings, ...
            'zero_axis_all_zero',o.zero_axis_all_zero, ...
            'amplitude_asymmetry_l2',NaN, ...
            'upper_background_median_tl_db',median(values), ...
            'upper_background_valid_fraction',nnz(isfinite(tl(background)))/nnz(background), ...
            'upper_background_fraction_ge70',nnz(isfinite(tl(background)) & tl(background)>=70)/nnz(background))]; %#ok<AGROW>
        % Relative amplitude symmetry is diagnostic only; no data are symmetrized.
        mirrored=fliplr(abs(o.pressure));
        summary(end).amplitude_asymmetry_l2=norm(abs(o.pressure(mask))-mirrored(mask))/norm(abs(o.pressure(mask)));
    end
    for pair=[1 2 3;2 3 4]
        a=fields{pair(1)}; b=fields{pair(2)};
        mask=a.depth_m>=10 & a.depth_m<=4900 & abs(a.range_m)>=.1;
        both=mask & abs(a.pressure)>0 & abs(b.pressure)>0;
        delta=sort(abs(-20*log10(abs(a.pressure(both)./b.pressure(both)))));
        if isempty(delta), q=[NaN NaN]; else, q=delta(ceil([.5 .95]*numel(delta))); end
        checks=[checks struct('attempt',attempt,'a',string(a.cfg.tag),'b',string(b.cfg.tag), ...
            'complex_relative_l2',norm(a.pressure(mask)-b.pressure(mask))/norm(b.pressure(mask)), ...
            'common_nonzero_fraction',nnz(both)/nnz(mask), ...
            'tl_median_absolute_db',q(1),'tl_p95_absolute_db',q(2))]; %#ok<AGROW>
    end
    single=figure('Visible','off','Color','w','Position',[40 40 780 840]);
    ax=axes('Parent',single,'Units','normalized','Position',[.14 .1 .69 .82]);
    local_field(ax,fields{4}); title(ax,labels{attempt},'FontSize',13,'Interpreter','none');
    ax.Position=[.14 .1 .69 .82]; cb=findobj(single,'Type','ColorBar');
    cb.Position=[.85 .1 .025 .82];
    local_export(single,fullfile(figdir,sprintf('attempt%d',attempt)));
    s=load(fullfile(outdir,sprintf('a%d_control.mat',attempt)),'output'); o=s.output;
    distance=sqrt((o.cfg.source_depth_m-o.depth_m).^2+o.range_m.^2);
    if attempt<3, reference=1./distance;
    else, reference=4*sqrt(pi)./sqrt(distance); end
    for dd=1:numel(o.depth_m)
        mask=abs(o.range_m)>=.1 & abs(o.range_m)<=5;
        nz=mask & abs(o.pressure(dd,:))>0;
        err=sort(abs(-20*log10(abs(o.pressure(dd,nz))./reference(dd,nz))));
        if isempty(err), q=[NaN NaN]; else, q=err(ceil([.5 .95]*numel(err))); end
        controls=[controls struct('attempt',attempt,'depth_m',o.depth_m(dd), ...
            'common_nonzero_fraction',nnz(nz)/nnz(mask), ...
            'tl_median_error_db',q(1),'tl_p95_error_db',q(2))]; %#ok<AGROW>
    end
end
title(layout,'Three bounded native attempts: coherent TL, 20-70 dB; gray = zero/invalid');
local_export(fig,fullfile(figdir,'three_attempt_comparison'));
summary=struct2table(summary); checks=struct2table(checks); controls=struct2table(controls);
writetable(summary,fullfile(outdir,'case_summary.csv'));
writetable(checks,fullfile(outdir,'convergence.csv'));
writetable(controls,fullfile(outdir,'uniform_control.csv'));
save(fullfile(outdir,'assessment.mat'),'summary','checks','controls','labels');
disp(summary); disp(checks); disp(controls);
end

function local_field(ax,o)
tl=-20*log10(abs(o.pressure)); valid=isfinite(tl);
% Explicit display RGB preserves gray gaps without renderer transparency.
map=flipud(jet(256)); index=1+floor(255*min(max((tl-20)/50,0),1));
index(~valid)=1; rgb=ind2rgb(index,map);
for component=1:3
    plane=rgb(:,:,component); plane(~valid)=.7; rgb(:,:,component)=plane;
end
image(ax,o.range_m,o.depth_m,rgb);
ax.Color=[.7 .7 .7]; ax.YDir='reverse'; xlim(ax,[-5 5]); ylim(ax,[0 5000]);
colormap(ax,map); clim(ax,[20 70]);
bar=colorbar(ax); bar.Direction='reverse'; bar.Label.String='TL (dB)';
xlabel(ax,'Range (m)','Interpreter','none'); ylabel(ax,'Depth (m)','Interpreter','none');
ax.FontSize=12; ax.FontName='Arial'; bar.FontName='Arial';
end

function local_export(fig,path)
drawnow;
% Export the whole canvas: avoid automatic crop cutting labels in this session.
fig.PaperUnits='inches'; size_inch=fig.Position(3:4)/get(groot,'ScreenPixelsPerInch');
fig.PaperPosition=[0 0 size_inch]; fig.PaperSize=size_inch;
print(fig,[path '.png'],'-dpng','-vector','-r180');
print(fig,[path '.pdf'],'-dpdf','-vector','-r180'); close(fig);
end
