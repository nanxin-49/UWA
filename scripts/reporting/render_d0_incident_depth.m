function receipt = render_d0_incident_depth(out,cfg,label)
%RENDER_D0_INCIDENT_DEPTH Saved-data figures only. Does not evaluate a source.
arguments
    out (1,:) char
    cfg (1,1) struct
    label (1,:) char = 'D0 INCIDENT ONLY'
end
folder=fullfile(out,'figures'); assert(~isfolder(folder)); mkdir(folder);
all_data=cell(5,1); metrics=cell(5,1);
for j=1:5
    v=load(fullfile(out,sprintf('depth_%04dm',cfg.depths_m(j)),'field.mat'),'data');
    all_data{j}=v.data; metrics{j}=v.data.metrics;
end
writetable(struct2table(vertcat(metrics{:})),fullfile(folder,'depth_summary.csv'));
colors=lines(5); styles={'-','--',':','-.','-'};
titles={'Absolute incident peak pressure','Two explicit normalizations', ...
    'Phase relative to physical center value','Finite-window and continuous angular spectra', ...
    'Geometric rays and local phase directions','Pressure-gradient magnitude and phase curvature'};
stems={'F01_absolute_pressure','F02_normalized_amplitude','F03_relative_phase', ...
    'F04_complex_angular_spectrum','F05_incidence_angles','F06_gradients_curvature'};
total_lines=0;
for fi=1:6
    fig=figure('Visible','off','Position',[60 60 1200 820],'Color','w');
    guard=onCleanup(@()close(fig));
    if fi==4, np=4; else, np=1+ismember(fi,[2 5 6]); end
    if np==4, tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
    else, tl=tiledlayout(fig,np,1,'Padding','compact','TileSpacing','compact'); end
    ax=gobjects(np,1);
    for panel=1:np
        ax(panel)=nexttile(tl); hold(ax(panel),'on'); grid(ax(panel),'on');
        for j=1:5
            d=all_data{j}; x=d.x_m;
            switch fi
                case 1, y=abs(d.p_f); yl='Peak pressure (Pa)';
                case 2
                    if panel==1, y=abs(d.center_normalized_f); yl='|p / p(0)|';
                    else, y=abs(d.rms_normalized_f); yl='|p| / RMS line pressure'; end
                case 3, y=d.relative_phase_rad; yl='Unwrapped phase - phase at x=0 (rad)';
                case 4
                    if panel<=2
                        w=d.window_spectra{4}; x=w.kx_rad_per_m;
                        if panel==1, y=abs(w.spectrum_f); yl='440 m Hann window |P| (Pa m)';
                        else, y=w.phase_rad; yl='440 m Hann window wrapped phase (rad)'; end
                    else
                        q=d.continuous_spectrum; x=q.prop_kx_rad_per_m;
                        if panel==3, y=abs(q.prop_density_f); yl='Continuous propagating |P| (Pa m)';
                        else, y=q.prop_phase_unwrapped_rad; yl='Continuous unwrapped phase (rad)'; end
                    end
                case 5
                    if panel==1, y=d.geometric_angle_deg; yl='Center-ray angle to upward vertical (deg)';
                    else, y=d.local_angle_deg; yl='Local phase-gradient direction (deg)'; end
                case 6
                    if panel==1, y=hypot(abs(d.px_f),abs(d.pz_f)); yl='Gradient norm (Pa/m)';
                    else, y=d.phase_curvature_rad_per_m2; yl='d(phase k_x)/dx (rad/m^2)'; end
            end
            line_handle=plot(ax(panel),x,y,'Color',colors(j,:),'LineStyle',styles{j}, ...
                'LineWidth',1.4,'DisplayName',sprintf('Tx %g m',cfg.depths_m(j)));
            if j==5, set(line_handle,'Marker','o','MarkerIndices',1:1000:numel(x),'MarkerSize',4); end
        end
        % Do not visually magnify floating-point noise in constant fixtures/fields.
        owned=findall(ax(panel),'Type','line'); values=[];
        for ll=1:numel(owned), values=[values owned(ll).YData]; end %#ok<AGROW>
        values=values(isfinite(values));
        if ~isempty(values)
            if fi==5
                ylim(ax(panel),[-90 90]);
            elseif fi==1 || fi==2 || (fi==6 && panel==1)
                ylim(ax(panel),[0 max(max(values)*1.08,1e-6)]);
            elseif max(values)-min(values)<1e-8*max(1,max(abs(values)))
                middle=mean(values); pad=max(1e-3,.02*abs(middle));
                ylim(ax(panel),middle+[-pad pad]);
            end
        end
        ylabel(ax(panel),yl,'Interpreter','none');
        if fi==4
            xlabel(ax(panel),'Horizontal wavenumber k_x (rad/m)');
            xlim(ax(panel),[-1 1]*(2*pi*all_data{1}.source.frequency_hz/all_data{1}.source.c_mps)*1.05);
            if panel==1 || panel==3, set(ax(panel),'YScale','log'); end
        else
            xlabel(ax(panel),'Horizontal position x (m)');
        end
    end
    title(tl,[label ' | ' titles{fi}],'Interpreter','none');
    lg=legend(ax(1),'Orientation','horizontal'); lg.Layout.Tile='south';
    drawnow; path=fullfile(folder,stems{fi}); savefig(fig,[path '.fig']);
    exportgraphics(fig,[path '.png'],'Resolution',160);
    % Reopen trusted FIG, compare every data line and labels, and reconstruct PNG.
    orig=findall(fig,'Type','line'); signature=local_signatures(orig);
    restored=openfig([path '.fig'],'new','invisible'); rg=onCleanup(@()close(restored));
    received=local_signatures(findall(restored,'Type','line'));
    assert(isequaln(signature,received),'FIG data changed');
    restored.Position=fig.Position; drawnow;
    exportgraphics(restored,[path '_reopened.png'],'Resolution',160);
    assert(all(isfile({[path '.fig'],[path '.png'],[path '_reopened.png']})));
    total_lines=total_lines+numel(orig); clear rg guard
end
receipt=struct('status','FIG_DATA_REOPEN_PASS','figures',6,'verified_line_count',total_lines, ...
    'failed_checks',0,'PNG_reconstruction','exported from reopened FIG; visual inspection after Fetch remains required', ...
    'pixel_equality_required',false,'label',label);
pm_sensitivity_io('write',fullfile(folder,'figure_audit.json'),receipt);
fid=fopen(fullfile(out,'D0_run_report.md'),'w','n','UTF-8'); assert(fid>=0); guard=onCleanup(@()fclose(fid));
fprintf(fid,'# %s: saved incident-field receipt\n\n',label);
fprintf(fid,'Five depths: 100/300/600/1000/1490 m. Width .10 m; 4 kHz; v0=.001 m/s.\n\n');
fprintf(fid,'Uniform lossless c=1500 m/s, rho=1025 kg/m^3; exp(-i omega t), absolute physical Pa, phase contract v2.\n\n');
fprintf(fid,'No bottom boundary, reflection, images, shell, downward radiation, SSP, absorption or current. This isolates source-depth effects, not a complete real deep-sea environment.\n\n');
fprintf(fid,'See figures/depth_summary.csv, depth_*/audit_checks.csv and plane_wave.csv. Finite-window FFT outside-band power is not physical evanescent energy. Continuous grazing singularities are omitted from plotted density, retained through transformed integration in independent audits.\n\n');
fprintf(fid,'K0 U_PM=6 m/s, N=12, seed260001, L160 m, original 42--50 m C2. Rough evaluations only at original saved coordinates inside +/-50 m. No tiling or extension; D1 coverage is unresolved.\n\n');
fprintf(fid,'Center normalization divides by complex p(0); RMS normalization divides by a positive whole-line RMS and retains physical phase. Relative-phase plots subtract only the declared center reference. No fitted phase. Tangent plane waves use actual anchor pressure and local phase direction; inspect their saved patch errors, not distance alone.\n\n');
fprintf(fid,'R1/R2 Tx100 m reflection errors cannot be transferred to Tx1490 m without new paired reflection references. This run contains incident fields only. Final physical interpretation and PNG visual review await manual Fetch and Codex analysis.\n');
end
function signatures=local_signatures(lines_in)
signatures=cell(numel(lines_in),1);
for j=1:numel(lines_in)
    h=lines_in(j); signatures{j}=struct('x',h.XData,'y',h.YData,'name',h.DisplayName, ...
        'xlabel',h.Parent.XLabel.String,'ylabel',h.Parent.YLabel.String);
end
% Traversal order is stable after savefig/openfig for this explicit object tree.
end
