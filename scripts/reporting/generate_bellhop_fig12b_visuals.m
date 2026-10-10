function files = generate_bellhop_fig12b_visuals(output_dir, scope)
%GENERATE_BELLHOP_FIG12B_VISUALS Plot native saved fields, retaining gaps.
arguments
    output_dir (1,:) char = ''
    scope (1,:) char = 'all'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(output_dir)
    output_dir=fullfile(root,'results','validation','bellhop_fig12b');
end
figure_dir=fullfile(output_dir,'figures');
if ~isfolder(figure_dir), mkdir(figure_dir); end
if strcmpi(scope,'rays')
    local_rayfig(output_dir,figure_dir);
    files=struct('figure_dir',figure_dir); return
end
assert(strcmpi(scope,'all'),'scope must be all or rays.');
listing=dir(fullfile(output_dir,'*.mat'));
cases={}; rows=struct([]);
for ii=1:numel(listing)
    if strcmp(listing(ii).name,'plan.mat'), continue; end
    s=load(fullfile(output_dir,listing(ii).name),'output'); o=s.output;
    if isempty(o.pressure), continue; end
    cases{end+1}=o; %#ok<AGROW>
    finite=isfinite(o.pressure); nz=finite & abs(o.pressure)>0;
    tl=-20*log10(abs(o.pressure(nz)));
    if isempty(tl), tl=NaN; end
    rows=[rows struct('tag',string(o.cfg.tag),'beams',o.cfg.beam_count, ...
        'step_m',o.cfg.step_m,'seconds',o.elapsed_s, ...
        'zero_fraction',o.zero_fraction,'nonfinite_count',o.nonfinite_count, ...
        'storage_warnings',o.storage_warnings,'unbounded_warnings',o.unbounded_warnings, ...
        'tl_min_db',min(tl),'tl_max_db',max(tl))]; %#ok<AGROW>
end
writetable(struct2table(rows),fullfile(output_dir,'case_summary.csv'));
files=struct('figure_dir',figure_dir);
for fan={'initial','bilateral','wide'}
    fig=figure('Visible','off','Color','w','Position',[50 50 1350 610]);
    layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    for n=[50 500 2000]
        tag=sprintf('%s_n%d_h10',fan{1},n);
        o=local_find(cases,tag);
        ax=nexttile(layout); local_field(ax,o);
        title(ax,sprintf('%d beams; step 10 m',n));
    end
    title(layout,sprintf('%s fan: native coherent point-source TL; gray = zero / invalid',fan{1}));
    local_export(fig,fullfile(figure_dir,[fan{1} '_beam_comparison']));
end
fig=figure('Visible','off','Color','w','Position',[50 50 1350 610]);
layout=tiledlayout(fig,1,4,'TileSpacing','compact','Padding','compact');
for h=[10 5 2.5 1.25]
    o=local_find(cases,sprintf('wide_n2000_h%g',h));
    ax=nexttile(layout); local_field(ax,o); title(ax,sprintf('step %g m',h));
end
title(layout,'Wide fan, 2000 beams: integration-step comparison');
local_export(fig,fullfile(figure_dir,'wide_step_comparison'));
comparisons={ ...
    'initial_n50_h10','initial_n500_h10'; ...
    'initial_n500_h10','initial_n2000_h10'; ...
    'bilateral_n50_h10','bilateral_n500_h10'; ...
    'bilateral_n500_h10','bilateral_n2000_h10'; ...
    'wide_n50_h10','wide_n500_h10'; ...
    'wide_n500_h10','wide_n2000_h10'; ...
    'wide_n2000_h10','wide_n2000_h5'; ...
    'wide_n2000_h5','wide_n2000_h2.5'; ...
    'wide_n2000_h2.5','wide_n2000_h1.25'; ...
    'wide_n2000_h1.25','wide_n2000_h0.625'; ...
    'wide_n2000_h0.625','wide_n2000_h0.3125'; ...
    'wide_n2000_h0.3125','wide_n10000_h0.3125'; ...
    'sym_n50_h0.3125','sym_n500_h0.3125'; ...
    'sym_n500_h0.3125','sym_n2000_h0.3125'; ...
    'sym_n2000_h0.3125','sym_n2000_h0.15625'; ...
    'bilateral_n2000_h10','bilateral_n2000_h1.25'; ...
    'bilateral_n2000_h1.25','bilateral_n2000_h0.625'; ...
    'bilateral_n2000_h0.625','bilateral_n2000_h0.3125'; ...
    'bilateral_n500_h0.625','bilateral_n2000_h0.625'; ...
    'bilateral_n2000_h0.625','bilateral_n8000_h0.625'};
metrics=struct([]);
for ii=1:size(comparisons,1)
    a=local_find(cases,comparisons{ii,1}); b=local_find(cases,comparisons{ii,2});
    % Fixed physical mask; no phase alignment, fitted gain or color clipping.
    mask=a.depth_m>=10 & a.depth_m<=4900 & abs(a.range_m)>=0.1;
    nonzero=abs(a.pressure)>0 & abs(b.pressure)>0;
    both=mask & nonzero & isfinite(a.pressure) & isfinite(b.pressure);
    delta=abs(-20*log10(abs(a.pressure(both)./b.pressure(both))));
    delta=sort(delta); q=delta(max(1,ceil([.5 .95]*numel(delta))));
    metrics=[metrics struct('a',string(comparisons{ii,1}), ...
        'b',string(comparisons{ii,2}), ...
        'complex_relative_l2',norm(a.pressure(mask)-b.pressure(mask))/norm(b.pressure(mask)), ...
        'common_nonzero_fraction',nnz(both)/nnz(mask), ...
        'tl_median_absolute_db',q(1),'tl_p95_absolute_db',q(2))]; %#ok<AGROW>
end
writetable(struct2table(metrics),fullfile(output_dir,'convergence.csv'));
disp(struct2table(metrics));
for tag={'bilateral_n2000_h0.625','sym_n2000_h0.15625','bilateral_n8000_h0.625'}
    o=local_find(cases,tag{1});
    fig=figure('Visible','off','Color','w','Position',[50 50 750 830]);
    ax=nexttile(tiledlayout(fig,1,1)); local_field(ax,o);
    title(ax,sprintf('%d beams, step %g m; [%g,%g] deg', ...
        o.cfg.beam_count,o.cfg.step_m,o.cfg.angles_deg));
    local_export(fig,fullfile(figure_dir,tag{1}));
end
fig=figure('Visible','off','Color','w','Position',[50 50 1350 610]);
layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for n=[50 500 2000]
    o=local_find(cases,sprintf('sym_n%d_h0.3125',n));
    ax=nexttile(layout); local_field(ax,o); title(ax,sprintf('%d symmetric beams',n));
end
title(layout,'Native full-angle comparison with paired left/right launch angles');
local_export(fig,fullfile(figure_dir,'symmetric_wide_beam_comparison'));
local_controls(cases,figure_dir,output_dir);
local_rayfig(output_dir,figure_dir);
end

function o=local_find(cases,tag)
for ii=1:numel(cases)
    if strcmp(cases{ii}.cfg.tag,tag), o=cases{ii}; return; end
end
error('Missing saved case %s',tag);
end

function local_controls(cases,figure_dir,output_dir)
fig=figure('Visible','off','Color','w','Position',[50 50 1300 620]);
layout=tiledlayout(fig,1,2); rows=struct([]);
for jj=1:2
    tags={'control_cart_2000','control_narrow_2000'};
    o=local_find(cases,tags{jj}); ax=nexttile(layout); hold(ax,'on');
    distance=sqrt((o.cfg.source_depth_m-o.depth_m).^2+o.range_m.^2);
    ratio=abs(o.pressure).*distance;
    for dd=1:numel(o.depth_m)
        plot(ax,o.range_m,ratio(dd,:),'DisplayName',sprintf('z=%g m',o.depth_m(dd)));
    end
    yline(ax,1,'k--','1/R reference','HandleVisibility','off');
    xlabel(ax,'Range (m)'); ylabel(ax,'|p| R (no fitted gain)');
    title(ax,tags{jj},'Interpreter','none'); legend(ax,'Location','best');
    xlim(ax,[-5 5]); grid(ax,'on');
    for dd=1:numel(o.depth_m)
        % Interior of the directly illuminated fan, away from its edges/axis.
        dz=o.cfg.source_depth_m-o.depth_m(dd);
        if jj==1, limit=5; else, limit=dz*tand(.055); end
        mask=abs(o.range_m)>=.1 & abs(o.range_m)<=limit;
        delta=abs(20*log10(ratio(dd,mask))); delta=sort(delta);
        if isempty(delta), med=NaN; p95=NaN;
        else, med=delta(ceil(.5*numel(delta))); p95=delta(ceil(.95*numel(delta))); end
        rows=[rows struct('tag',string(tags{jj}),'depth_m',o.depth_m(dd), ...
            'samples',nnz(mask),'median_tl_error_db',med,'p95_tl_error_db',p95)]; %#ok<AGROW>
    end
end
title(layout,'Independent matched-halfspace homogeneous control: exact spherical magnitude 1/R');
local_export(fig,fullfile(figure_dir,'analytic_control'));
writetable(struct2table(rows),fullfile(output_dir,'analytic_control.csv'));
o=local_find(cases,'surface_patch');
fig=figure('Visible','off','Color','w','Position',[50 50 900 570]);
ax=axes(fig); imagesc(ax,o.range_m,o.depth_m,-20*log10(abs(o.pressure)), ...
    'AlphaData',abs(o.pressure)>0); set(ax,'YDir','reverse','Color',[.75 .75 .75]);
colormap(ax,flipud(jet(256))); clim(ax,[20 90]); cb=colorbar(ax); cb.Direction='reverse';
xlabel(ax,'Range (m)'); ylabel(ax,'Depth (m)');
title(ax,'Near-surface coherent samples: dz=0.01 m (full-depth figures use 5 m)');
local_export(fig,fullfile(figure_dir,'surface_patch'));
end

function local_rayfig(output_dir,figure_dir)
fig=figure('Visible','off','Color','w','Position',[50 50 1300 680]);
layout=tiledlayout(fig,1,2); rows=struct([]);
for jj=1:2
    tags={'rays_initial','rays_bilateral'};
    s=load(fullfile(output_dir,[tags{jj} '.mat']),'output'); o=s.output;
    ax=nexttile(layout); hold(ax,'on');
    for kk=1:numel(o.rays)
        ray=o.rays(kk); rz=ray.rz_m;
        % Show up to the first return to the source depth after a surface hit.
        itop=find(rz(:,2)<1e-6,1);
        last=size(rz,1);
        if ~isempty(itop)
            ret=find(rz(itop:end,2)>=4950,1);
            if ~isempty(ret), last=itop+ret-1; end
        end
        selected=unique(round(linspace(1,numel(o.rays),12)));
        if ismember(kk,selected)
            % Display decimation avoids sub-pixel duplicate-segment rendering.
            % Keep the surface turning point and neighbors; raw RAY is intact.
            display_ix=unique([1:20:last last]);
            if ~isempty(itop)
                neighbors=max(1,itop-1):min(last,itop+1);
                display_ix=unique([display_ix neighbors]);
            end
            plot(ax,rz(display_ix,1),rz(display_ix,2),'Color',[.15 .4 .7], ...
                'LineStyle','-','LineWidth',.8,'Marker','none');
        end
        rows=[rows struct('tag',string(tags{jj}),'angle_deg',ray.angle_deg, ...
            'top_bounces',ray.top_bounces,'bottom_bounces',ray.bottom_bounces, ...
            'min_z_m',min(rz(:,2)),'max_z_m',max(rz(:,2)), ...
            'min_range_m',min(rz(:,1)),'max_range_m',max(rz(:,1)))]; %#ok<AGROW>
    end
    plot(ax,0,4950,'ro','MarkerFaceColor','r');
    set(ax,'YDir','reverse'); xlim(ax,[-5 5]); ylim(ax,[0 5000]);
    xlabel(ax,'Range (m)'); ylabel(ax,'Depth (m)'); grid(ax,'on');
    title(ax,tags{jj},'Interpreter','none');
end
title(layout,'12 representative trajectories from each native 50-ray run; no mirroring');
local_export(fig,fullfile(figure_dir,'representative_rays'));
writetable(struct2table(rows),fullfile(output_dir,'ray_audit.csv'));
end
function local_field(ax,o)
tl=-20*log10(abs(o.pressure));
valid=isfinite(tl);
imagesc(ax,o.range_m,o.depth_m,tl,'AlphaData',valid);
set(ax,'YDir','reverse','Color',[.75 .75 .75],'FontSize',12);
colormap(ax,flipud(jet(256))); clim(ax,[20 70]);
xlim(ax,[-5 5]); ylim(ax,[0 5000]);
xlabel(ax,'Range (m)'); ylabel(ax,'Depth (m)');
cb=colorbar(ax); cb.Direction='reverse'; cb.Label.String='TL (dB re native source)';
end

function local_export(fig,stem)
exportgraphics(fig,[stem '.png'],'Resolution',180);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
close(fig);
end
