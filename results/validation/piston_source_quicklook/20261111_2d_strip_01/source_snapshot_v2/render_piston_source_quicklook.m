function receipt = render_piston_source_quicklook(run_dir,figure_dir)
%RENDER_PISTON_SOURCE_QUICKLOOK Plot trusted saved 2-D aperture pressures.
% Reconstructs every figure from MAT data; no acoustic solver is invoked.
arguments
    run_dir (1,:) char
    figure_dir (1,:) char = fullfile(run_dir,'figures_final')
end
assert(isfolder(run_dir));
assert(~isfolder(figure_dir),'Use a fresh figure directory.'); mkdir(figure_dir);
colors=[.00 .35 .65;.80 .32 .00;.35 .35 .35]; styles={'-','--',':'};
labels={'20 cm total width','10 cm total width (historical 5 cm half-size)','5 cm total width'};
map_stems={};
for fi=1:3
    fkhz=[4 6 8]; s=load(fullfile(run_dir,sprintf('field_%dkHz.mat',fkhz(fi))));
    for other=[3 2]
        fig=figure('Visible','off','Color','w','Position',[40 40 1420 880]);
        tl=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('%d kHz | Exact 2-D baffled strip | v_0=1 mm/s (peak) | widths are TOTAL',fkhz(fi)),'FontSize',16);
        subtitle(tl,'Absolute levels: common [0,70] dB re 1 Pa. Ratios: common [-24,+24] dB. Display limits only; raw data retained.','FontSize',10);
        for row=1:2
            if row==1
                xx=s.x_m; zz=s.z_m; pp=s.p_field_f;
                ylimits=[0 100]; xlimits=[-120 120];
            else
                xx=s.x_near_m; zz=100-s.s_near_m; pp=s.p_near_f;
                ylimits=[99 100]; xlimits=[-.4 .4];
            end
            for col=1:3
                if row==1, ax=nexttile(tl,col); else, ax=nexttile(tl,col+3,[2 1]); end
                hold(ax,'on');
                if col<3
                    widx=[1 other]; wi=widx(col);
                    imagesc(ax,xx,zz,20*log10(abs(pp(:,:,wi))));
                    clim(ax,[0 70]); colormap(ax,parula(256));
                    title(ax,sprintf('W = %g cm | absolute',100*s.cfg.total_widths_m(wi)));
                    cb=colorbar(ax); cb.Label.String='Peak level (dB re 1 Pa)';
                else
                    imagesc(ax,xx,zz,20*log10(abs(pp(:,:,1))./abs(pp(:,:,other))));
                    clim(ax,[-24 24]); colormap(ax,local_diverge());
                    title(ax,sprintf('20 log_{10}(|p_{20}| / |p_{%g}|)',100*s.cfg.total_widths_m(other)));
                    cb=colorbar(ax); cb.Label.String='Amplitude ratio (dB)';
                end
                set(ax,'YDir','reverse','FontSize',11,'Layer','top');
                axis(ax,'image'); xlim(ax,xlimits); ylim(ax,ylimits);
                xlabel(ax,'Horizontal x (m)'); ylabel(ax,'Depth z (m), positive down');
                if row==1
                    h=yline(ax,3,'w--','z=3 m'); h.HandleVisibility='off';
                    plot(ax,0,100,'k^','MarkerFaceColor','w');
                else
                    text(ax,-.38,99.96,'Source plane z=100 m; s < 1 cm omitted', ...
                        'FontSize',9,'BackgroundColor','w','VerticalAlignment','top');
                end
            end
        end
        stem=sprintf('field_%dkHz_20vs%02dcm',fkhz(fi),round(100*s.cfg.total_widths_m(other)));
        local_save(fig,figure_dir,stem); map_stems{end+1}=stem; %#ok<AGROW>
    end
end
% One column per frequency: absolute Pa, normalized reception, aperture ratio.
fig=figure('Visible','off','Color','w','Position',[40 40 1440 970]);
tl=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
title(tl,'Receiver line z=3 m | 97 m above the source | exact aperture pressure','FontSize',17);
subtitle(tl,'Row 1: absolute peak Pa. Row 2: each curve divided by its own x=0 value. Gray lines bound |x|<=40 m.','FontSize',11);
for fi=1:3
    fkhz=[4 6 8]; s=load(fullfile(run_dir,sprintf('field_%dkHz.mat',fkhz(fi))));
    for row=1:3
        ax=nexttile(tl,(row-1)*3+fi); hold(ax,'on');
        for wi=1:3
            if row==1
                yy=abs(s.p_rx_f(:,wi));
            elseif row==2
                yy=20*log10(abs(s.p_rx_f(:,wi))/abs(s.p_rx_f(s.x_rx_m==0,wi)));
            else
                if wi==1, continue; end
                yy=20*log10(abs(s.p_rx_f(:,1))./abs(s.p_rx_f(:,wi)));
            end
            plot(ax,s.x_rx_m,yy,styles{wi},'Color',colors(wi,:),'LineWidth',1.9,'DisplayName',labels{wi});
        end
        grid(ax,'on'); xlim(ax,[-120 120]);
        h=xline(ax,[-40 40],'-','Color',[.65 .65 .65]); set(h,'HandleVisibility','off');
        xlabel(ax,'Horizontal x (m)');
        if row==1
            ylim(ax,[0 105]); ylabel(ax,'|p| (peak Pa)'); title(ax,sprintf('%d kHz',fkhz(fi)));
        elseif row==2
            ylim(ax,[-18 1]); ylabel(ax,'Axis-normalized amplitude (dB)');
        else
            ylim(ax,[-6 13]); ylabel(ax,'20 log_{10}(|p_{20}|/|p_W|) (dB)');
        end
        if fi==1 && row==1, legend(ax,'Location','south','FontSize',9); end
        if fi==1 && row==3
            legend(ax,{'20 / 10 cm','20 / 5 cm'},'Location','south','FontSize',10);
        end
        set(ax,'FontSize',11);
    end
end
local_save(fig,figure_dir,'receiver_line_comparison');
% Equal-radius circles separate directionality from cylindrical spreading.
fig=figure('Visible','off','Color','w','Position',[40 40 1440 760]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(tl,'2-D strip directivity | same radius 97 m | axis-normalized amplitude','FontSize',17);
subtitle(tl,'Upper: exact Hankel aperture integral (lines), far-field sinc (sparse markers). Lower: near-axis view. Normalized plots do not compare absolute output.','FontSize',11);
for fi=1:3
    fkhz=[4 6 8]; s=load(fullfile(run_dir,sprintf('field_%dkHz.mat',fkhz(fi))));
    for row=1:2
        ax=nexttile(tl,(row-1)*3+fi); hold(ax,'on');
        for wi=1:3
            dn=abs(s.p_circle_f(:,wi))/abs(s.p_circle_f(s.theta_deg==0,wi));
            plot(ax,s.theta_deg,20*log10(dn),styles{wi},'Color',colors(wi,:),'LineWidth',1.9,'DisplayName',labels{wi});
            ii=1:24:numel(s.theta_deg);
            plot(ax,s.theta_deg(ii),20*log10(abs(s.D_sinc(ii,wi))), ...
                'o','Color',colors(wi,:),'MarkerSize',3,'HandleVisibility','off');
        end
        h=yline(ax,-3,'k:','-3 dB'); h.HandleVisibility='off';
        h=xline(ax,atan2d(40,97)*[-1 1],'-','Color',[.7 .7 .7]); set(h,'HandleVisibility','off');
        if row==1, xlim(ax,[-90 90]); ylim(ax,[-45 1]); title(ax,sprintf('%d kHz',fkhz(fi)));
        else, xlim(ax,[-25 25]); ylim(ax,[-4 0.2]); end
        grid(ax,'on'); xlabel(ax,'Angle from upward vertical (deg)'); ylabel(ax,'20 log_{10}(|p|/|p_{axis}|) (dB)');
        set(ax,'FontSize',11);
        if fi==1 && row==1, legend(ax,'Location','south','FontSize',9); end
    end
end
local_save(fig,figure_dir,'directivity_comparison');
receipt=struct('status','FIGURES_RECONSTRUCTED_FROM_MAT','png_count',8,'fig_count',8,'field_stems',{map_stems},'output_dir',figure_dir,'data_dir',run_dir);
fid=fopen(fullfile(figure_dir,'render_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt));
end
function local_save(fig,out,stem)
assert(~isfile(fullfile(out,[stem '.png'])) && ~isfile(fullfile(out,[stem '.fig'])),'Refusing to overwrite figures.');
drawnow; savefig(fig,fullfile(out,[stem '.fig']));
exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',150); close(fig);
end
function c=local_diverge()
c=interp1([0 .5 1],[.12 .32 .65;1 1 1;.7 .12 .08],linspace(0,1,257));
end
