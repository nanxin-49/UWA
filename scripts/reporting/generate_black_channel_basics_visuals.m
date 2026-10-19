function output = generate_black_channel_basics_visuals(analysis,output_dir)
%GENERATE_BLACK_CHANNEL_BASICS_VISUALS Render saved BLACK analysis without rerun.
if ~isfolder(output_dir), mkdir(output_dir); end
representatives = analysis.cfg.representatives;
colors = [0 0.447 0.741;0.85 0.325 0.098;0.35 0.55 0.15];
styles = {'-','--',':'};
names = {'black_01_time_delay','black_02_power_delay','black_03_doppler'};
output = struct('png',{{}},'pdf',{{}},'figure_count',3);
fig = figure('Visible','off','Color','w','Position',[50 50 1500 650]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
reference = max(abs(analysis.h_reconstructed(:)));
for ii = 1:3
    ax = nexttile(layout);
    magnitude = abs(squeeze(analysis.h_reconstructed(:,ii,:)));
    db = 20*log10(max(magnitude/reference,realmin));
    imagesc(ax,analysis.reconstructed_delay_s*1e3,analysis.reconstructed_time_s,db.');
    axis(ax,'xy'); clim(ax,[-30 0]); colormap(ax,parula(256));
    xlabel(ax,'Relative delay (ms)'); ylabel(ax,'Record time (s)');
    title(ax,sprintf('Receiver %d',representatives(ii)));
    ax.FontSize = 12;
end
c = colorbar(ax); c.Label.String='Magnitude (dB / shared peak)';
title(layout,'BLACK: reconstructed time-varying multipath (native snapshot rate)');
local_export(fig,output_dir,names{1});
clear cleanup
fig = figure('Visible','off','Color','w','Position',[50 50 1500 760]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
for mm = 1:8
    ax = nexttile(layout);
    db = 10*log10(max(analysis.pdp(:,mm)/max(analysis.pdp(:,mm)),realmin));
    plot(ax,analysis.delay_s*1e3,db,'Color',[0 0.35 0.65],'LineWidth',1.3);
    hold(ax,'on'); yline(ax,-25,'--','-25 dB','Color',[0.5 0.25 0.1]);
    ylim(ax,[-45 1]); xlim(ax,[min(analysis.delay_s) max(analysis.delay_s)]*1e3);
    grid(ax,'on'); xlabel(ax,'Relative delay (ms)'); ylabel(ax,'Power (dB / element peak)');
    title(ax,sprintf('Rx %d; RMS = %.3f ms',mm,analysis.metrics.rms_delay_s(mm)*1e3));
    ax.FontSize = 10;
end
title(layout,'BLACK: mean stored-tap power-delay profiles; moments use -25 dB support');
local_export(fig,output_dir,names{2});
clear cleanup
fig = figure('Visible','off','Color','w','Position',[50 50 1400 900]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout,[1 2]); hold(ax,'on');
labels = arrayfun(@(m) sprintf('Rx %d',m),representatives,'UniformOutput',false);
if analysis.phase.available
    for ii = 1:3
        plot(ax,analysis.phase.time_s,analysis.phase.frequency_hz(:,representatives(ii)), ...
            styles{ii},'Color',colors(ii,:),'LineWidth',1.3);
    end
    legend(ax,labels,'Location','best');
else
    text(ax,0.1,0.5,'Phase fields absent: common frequency cannot be recovered.');
end
grid(ax,'on'); xlabel(ax,'Record time (s)'); ylabel(ax,'Common frequency shift (Hz)');
title(ax,'Continuous phase slope: non-overlapping 0.1 s fits; bulk resampling included');
ax.FontSize = 12;
for view = 1:2
    ax = nexttile(layout); hold(ax,'on');
    for ii = 1:3
        power = analysis.residual_psd(:,representatives(ii));
        db = 10*log10(max(power/max(power),realmin));
        plot(ax,analysis.residual_frequency_hz,db,styles{ii},'Color',colors(ii,:),'LineWidth',1.3);
    end
    grid(ax,'on'); ylim(ax,[-80 3]); xlabel(ax,'Residual Doppler frequency (Hz)');
    ylabel(ax,'PSD (dB / element peak)'); legend(ax,labels,'Location','best');
    if view==1
        xlim(ax,[-1 1]*analysis.params.fs_time/2);
        title(ax,sprintf('Stored main tap; native Nyquist +/- %.1f Hz',analysis.params.fs_time/2));
    else
        xlim(ax,[-1 1]*min(5,analysis.params.fs_time/2));
        title(ax,'Same spectrum: central +/- 5 Hz detail');
    end
    ax.FontSize = 12;
end
title(layout,'BLACK Doppler: tracked common frequency and residual tap spectrum are distinct');
local_export(fig,output_dir,names{3});
for ii=1:3
    output.png{ii} = fullfile(output_dir,[names{ii} '.png']);
    output.pdf{ii} = fullfile(output_dir,[names{ii} '.pdf']);
end
end

function local_export(fig,output_dir,name)
png = fullfile(output_dir,[name '.png']);
pdf = fullfile(output_dir,[name '.pdf']);
assert(~isfile(png) && ~isfile(pdf),'black:ExistingFigure','Refusing to overwrite figure.');
exportgraphics(fig,png,'Resolution',160);
exportgraphics(fig,pdf,'ContentType','vector');
end
