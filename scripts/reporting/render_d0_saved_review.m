function receipt = render_d0_saved_review(saved,out)
%RENDER_D0_SAVED_REVIEW Stable-layout scientific figures from saved D0 arrays.
folder=fullfile(out,'figures'); assert(~isfolder(folder)); mkdir(folder);
colors=[0 .447 .741;.85 .325 .098;.466 .674 .188;.494 .184 .556;0 .55 .55];
styles={'-','--',':','-.','-'}; labels=cellfun(@(d)sprintf('Tx %g m',d.source.z_source_m),saved,UniformOutput=false);
stems={'F01_absolute_pressure','F02_normalized_amplitudes','F03_relative_phase', ...
    'F04_angular_spectra','F05_incidence_angles','F06_gradients_curvature', ...
    'F07_tangent_plane_errors','F08_window_spectrum_phase'};
line_count=0;
for fi=1:8
    if fi==8
        hgt=1050; np=5;
    elseif fi==4, hgt=560; np=3;
    elseif ismember(fi,[2 5 6 7]), hgt=750; np=2;
    else
        hgt=570; np=1;
    end
    fig=figure('Visible','off','Color','w','Position',[40 40 1300 hgt]); fg=onCleanup(@()close(fig));
    axes_h=gobjects(np,1);
    for panel=1:np
        if fi==4
            pos=[.065+(panel-1)*.325 .23 .265 .64];
        elseif np==1, pos=[.10 .22 .86 .66];
        else
            gap=.055; total=.77; ph=(total-(np-1)*gap)/np; pos=[.11 .16+(np-panel)*(ph+gap) .85 ph];
        end
        ax=axes(fig,'Position',pos,'FontSize',11,'Box','on'); axes_h(panel)=ax; hold(ax,'on'); grid(ax,'on');
        for j=1:5
            d=saved{j}; x=d.x_m; q=d.continuous_spectrum;
            switch fi
                case 1, y=abs(d.p_f); yl='Peak pressure (Pa)';
                case 2
                    if panel==1, y=abs(d.center_normalized_f); yl='|p(x) / p(0)|';
                    else, y=abs(d.rms_normalized_f); yl='|p(x)| / spatial RMS'; end
                case 3, y=d.relative_phase_rad; yl='Unwrapped phase relative to x = 0 (rad)';
                case 4
                    if panel==1
                        w=d.window_spectra{4}; x=w.kx_rad_per_m; y=abs(w.spectrum_f); yl='440 m Hann FFT magnitude (Pa m)';
                    elseif panel==2, x=q.prop_kx_rad_per_m; y=abs(q.prop_density_f); yl='Continuous propagating density (Pa m)';
                    else
                        x=q.prop_kx_rad_per_m; y=q.prop_phase_unwrapped_rad; yl='Continuous physical phase (rad)';
                    end
                case 5
                    if panel==1, y=d.geometric_angle_deg; yl='Geometric incidence (deg)';
                    else, y=d.local_angle_deg; yl='Local phase direction (deg)'; end
                case 6
                    if panel==1, y=hypot(abs(d.px_f),abs(d.pz_f)); yl='Pressure gradient magnitude (Pa/m)';
                    else, y=d.phase_curvature_rad_per_m2; yl='Phase curvature (rad/m^2)'; end
                case 7
                    anchor=(panel-1)*200; t=d.plane_wave; sel=t.anchor_x_m==anchor;
                    x=t.patch_width_m(sel); y=100*t.complex_relative_L2(sel); yl='Complex relative L2 (%)';
                case 8
                    if j~=panel, continue; end
                    w=d.window_spectra{4}; x=w.kx_rad_per_m; y=w.phase_rad;
                    y(abs(x)>=2*pi*4000/1500)=NaN; yl=sprintf('Tx %g m\nPhase (rad)',d.source.z_source_m);
            end
            h=plot(ax,x,y,'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',1.5,'DisplayName',labels{j});
            if fi==7
                set(h,'Marker','o','MarkerSize',5);
            elseif fi==8, set(h,'LineStyle','none','Marker','.','MarkerSize',3);
            elseif j==5
                set(h,'Marker','o','MarkerIndices',1:1200:numel(x),'MarkerSize',4);
            end
        end
        ylabel(ax,yl,'Interpreter','none');
        if fi==4 || fi==8
            xlabel(ax,'Horizontal wavenumber k_x (rad/m)'); xlim(ax,[-17.5 17.5]);
            if fi==4 && panel<=2, set(ax,'YScale','log'); end
            if fi==8, ylim(ax,[-pi pi]); yticks(ax,[-pi 0 pi]); yticklabels(ax,{'-pi','0','pi'}); end
        elseif fi==7
            xlabel(ax,'Patch width (m)'); set(ax,'XScale','log','YScale','log'); xticks(ax,[1 5 10 40]); xlim(ax,[.9 45]);
            text(ax,.025,.90,sprintf('Anchor x_0 = %g m',(panel-1)*200),'Units','normalized', ...
                'FontSize',11,'BackgroundColor','w');
        else
            xlabel(ax,'Horizontal position x (m)'); xlim(ax,[-220 220]);
            if fi==1 || fi==2 || (fi==6 && panel==1), lim=ylim(ax); ylim(ax,[0 lim(2)]); end
        end
        if np>1 && fi~=4 && panel<np, ax.XLabel.String=''; end
    end
    titles={'D0: absolute incident pressure, 4 kHz / 0.10 m strip', ...
        'D0: center normalization and whole-line spatial RMS normalization', ...
        'D0: relative phase; no fitted depth-dependent phase offset', ...
        'D0: window spectrum versus full-line density (grazing endpoints excluded)', ...
        'D0: geometric rays and finite-aperture phase-gradient directions', ...
        'D0: analytic pressure gradients and phase curvature', ...
        'D0: tangent plane anchored to actual pressure and local phase direction', ...
        'D0: 440 m Hann FFT phase; valid stored bins, points without connecting wraps'};
    annotation(fig,'textbox',[.04 .935 .92 .055],'String',titles{fi},'EdgeColor','none', ...
        'HorizontalAlignment','center','FontSize',13,'FontWeight','bold','Interpreter','none');
    if fi~=8
        legend(axes_h(1),'Position',[.19 .045 .65 .045],'Orientation','horizontal','FontSize',10,'Box','off');
    end
    if fi==4
        text(axes_h(2),0,.9,'Five depth curves coincide in magnitude','Units','normalized','Interpreter','none', ...
            'HorizontalAlignment','left','FontSize',9);
    end
    drawnow; p=fullfile(folder,stems{fi}); savefig(fig,[p '.fig']); exportgraphics(fig,[p '.png'],'Resolution',160);
    original=signatures(fig); restored=openfig([p '.fig'],'new','invisible'); rg=onCleanup(@()close(restored));
    assert(isequaln(original,signatures(restored)),'Review FIG data/labels mismatch');
    restored.Position=fig.Position; drawnow; exportgraphics(restored,[p '_reopened.png'],'Resolution',160);
    line_count=line_count+numel(original); clear rg fg
end
receipt=struct('status','NEW_FIG_DATA_REOPEN_PASS','figures',8,'lines',line_count,'failed_checks',0, ...
    'input','Saved D0 arrays only','PNG_reconstruction','All eight re-exported from reopened trusted FIG; no cross-platform pixel gate');
pm_sensitivity_io('write',fullfile(folder,'figure_audit.json'),receipt);
end
function s=signatures(fig)
lines=findall(fig,'Type','line'); s=cell(numel(lines),1);
for j=1:numel(lines)
    s{j}=struct('x',lines(j).XData,'y',lines(j).YData,'name',lines(j).DisplayName, ...
        'xlabel',lines(j).Parent.XLabel.String,'ylabel',lines(j).Parent.YLabel.String);
end
end
