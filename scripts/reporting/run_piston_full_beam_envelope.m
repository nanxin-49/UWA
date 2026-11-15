function receipt = run_piston_full_beam_envelope(run_id)
%RUN_PISTON_FULL_BEAM_ENVELOPE Exact 2-D free-field pressure and ray sketches.
% Geometry rays are illustrative; pressure comes from the velocity aperture.
arguments
    run_id (1,:) char
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')),'Invalid run ID.');
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','visualization','piston_source_full_beam_envelope',run_id);
assert(~isfolder(out),'Result directory exists.'); mkdir(out);
logfile=fullfile(out,'matlab_diary.log'); diary(logfile); cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
t0=tic;
cfg=struct('frequencies_hz',[4000 6000 8000], ...
    'total_widths_m',[.20 .10],'source_x_m',0,'source_z_m',100, ...
    'receiver_z_m',3,'c_mps',1500,'rho_kgpm3',1025, ...
    'v0_peak_mps',.001,'field_nquad',24,'refined_nquad',96, ...
    'angle_radius_m',97,'angle_step_deg',.25, ...
    'pressure_units','peak Pa; exp(-i omega t)', ...
    'geometry','2-D infinite y strip in rigid baffle, free field', ...
    'ray_meaning','straight geometric directions, not wave paths or energy bounds');
x_m=-300:1:300; z_m=0:.5:99.5; [X,Z]=meshgrid(x_m,z_m);
theta_deg=(-89.5:cfg.angle_step_deg:89.5).';
p_field_f=complex(zeros(numel(z_m),numel(x_m),2,3));
p_circle_f=complex(zeros(numel(theta_deg),2,3));
D_sinc=zeros(numel(theta_deg),2,3);
quad_error=zeros(2,3); farfield_error=zeros(2,3);
half_3db_deg=nan(2,3); first_null_deg=nan(2,3);
ix_sample=round(linspace(1,numel(x_m),19));
iz_sample=unique([1 2 4 10 25 50 100 150 200]);
fprintf('RUN_ID=%s VERSION=%s ROOT=%s\n',run_id,version,root);
fprintf('GRID=%dx%d X=[%g,%g] Z=[%g,%g] NQUAD=%d\n', ...
    numel(z_m),numel(x_m),x_m(1),x_m(end),z_m(1),z_m(end),cfg.field_nquad);
try
    for fi=1:3
        f=cfg.frequencies_hz(fi); tf=tic;
        for wi=1:2
            W=cfg.total_widths_m(wi);
            source=struct('frequency_hz',f,'c_mps',cfg.c_mps, ...
                'rho_kgpm3',cfg.rho_kgpm3,'v0_mps',cfg.v0_peak_mps, ...
                'total_width_m',W,'z_source_m',cfg.source_z_m);
            p_field_f(:,:,wi,fi)=baffled_strip_piston_field_vertical(X,Z,source,cfg.field_nquad);
            [XS,ZS]=meshgrid(x_m(ix_sample),z_m(iz_sample));
            coarse=p_field_f(iz_sample,ix_sample,wi,fi);
            refined=baffled_strip_piston_field_vertical(XS,ZS,source,cfg.refined_nquad);
            quad_error(wi,fi)=max(abs(coarse(:)-refined(:)))/max(abs(refined(:)));
            assert(quad_error(wi,fi)<1e-7,'Quadrature refinement failed.');
            xc=cfg.angle_radius_m*sind(theta_deg);
            zc=cfg.source_z_m-cfg.angle_radius_m*cosd(theta_deg);
            p_circle_f(:,wi,fi)=baffled_strip_piston_field_vertical(xc,zc,source,cfg.refined_nquad);
            q=pi*f*W/cfg.c_mps*sind(theta_deg);
            d_sinc=ones(size(q)); nz=q~=0;
            d_sinc(nz)=sin(q(nz))./q(nz);
            D_sinc(:,wi,fi)=d_sinc;
            axis_idx=find(theta_deg==0,1);
            d=abs(p_circle_f(:,wi,fi))/abs(p_circle_f(axis_idx,wi,fi));
            farfield_error(wi,fi)=max(abs(d-abs(D_sinc(:,wi,fi))));
            assert(farfield_error(wi,fi)<1e-3,'Far-field sinc check failed.');
            positive=find(theta_deg>=0);
            j=find(d(positive)<=10^(-3/20),1);
            if ~isempty(j)
                k=positive(j); threshold=10^(-3/20);
                half_3db_deg(wi,fi)=theta_deg(k-1)+ ...
                    (theta_deg(k)-theta_deg(k-1))*(threshold-d(k-1))/(d(k)-d(k-1));
            end
            if f*W/cfg.c_mps>=1
                first_null_deg(wi,fi)=asind(cfg.c_mps/(f*W));
            end
            fprintf('F=%d W=%.2f QUAD_ERR=%.3g SINC_ERR=%.3g HALF_3DB=%.4g NULL=%.4g\n', ...
                f,W,quad_error(wi,fi),farfield_error(wi,fi),half_3db_deg(wi,fi),first_null_deg(wi,fi));
        end
        fprintf('F=%d FIELD_ELAPSED_S=%.3f\n',f,toc(tf));
    end
    save(fullfile(out,'complex_pressure_fields.mat'),'cfg','x_m','z_m','theta_deg', ...
        'p_field_f','p_circle_f','D_sinc','quad_error','farfield_error', ...
        'half_3db_deg','first_null_deg','-v7');
    clear X Z;
    cmax=ceil(max(20*log10(abs(p_field_f(:))))/10)*10;
    clim_pa=[0 cmax];
    for fi=1:3
        f=cfg.frequencies_hz(fi); fig=figure('Visible','off','Color','w','Position',[30 30 1500 1000]);
        tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('%d kHz | exact 2-D strip pressure above; geometric rays below | v_0=1 mm/s peak',f/1000),'FontSize',16);
        for wi=1:2
            W=cfg.total_widths_m(wi);
            ax=nexttile(tl,wi); hold(ax,'on');
            imagesc(ax,x_m,z_m,20*log10(max(abs(p_field_f(:,:,wi,fi)),realmin)));
            colormap(ax,parula(256)); clim(ax,clim_pa);
            cb=colorbar(ax); cb.Label.String='Absolute peak pressure (dB re 1 Pa)';
            plot(ax,0,100,'wv','MarkerFaceColor','k','MarkerSize',8,'HandleVisibility','off');
            yline(ax,3,'w--','z = 3 m','LabelHorizontalAlignment','left','HandleVisibility','off');
            title(ax,sprintf('Total strip width = %d cm | wave pressure',round(W*100)));
            format_axes(ax);
            ax=nexttile(tl,wi+2); hold(ax,'on');
            for th=-85:10:85
                [xe,ze,clipped]=ray_endpoint(th);
                if clipped, sty=':'; else, sty='-'; end
                plot(ax,[0 xe],[100 ze],sty,'Color',[.65 .69 .73],'LineWidth',.9,'HandleVisibility','off');
            end
            if isfinite(half_3db_deg(wi,fi))
                for sg=[-1 1]
                    th=sg*half_3db_deg(wi,fi); [xe,ze,clipped]=ray_endpoint(th);
                    plot(ax,[0 xe],[100 ze],'b--','LineWidth',2.2,'HandleVisibility','off');
                    if clipped, plot(ax,xe,ze,'bs','MarkerFaceColor','b','HandleVisibility','off'); end
                end
                note=sprintf('Dashed: +/-%.1f deg at -3 dB (97 m circle)',half_3db_deg(wi,fi));
            else
                note='No forward -3 dB boundary at 97 m; none drawn';
            end
            yline(ax,3,'k--','z = 3 m','LabelHorizontalAlignment','left','HandleVisibility','off');
            plot(ax,0,100,'kv','MarkerFaceColor',[.9 .35 .1],'MarkerSize',9,'HandleVisibility','off');
            text(ax,-288,93,sprintf('W=%d cm; actual source aperture +/-%.3f m',round(W*100),W/2), ...
                'FontSize',9,'BackgroundColor','w');
            text(ax,-288,83,note,'FontSize',9,'BackgroundColor','w');
            text(ax,175,20,'upward radiation','FontSize',10,'BackgroundColor','w');
            title(ax,sprintf('Total strip width = %d cm | ray fan schematic',round(W*100)));
            format_axes(ax);
        end
        stem=sprintf('full_beam_%dkHz',f/1000);
        savefig(fig,fullfile(out,[stem '.fig']));
        exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',180);
        close(fig);
    end
    fig=figure('Visible','off','Color','w','Position',[30 30 1500 550]);
    tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    title(tl,'Equal-radius exact 2-D directivity | 97 m from source | full forward angles','FontSize',16);
    for fi=1:3
        ax=nexttile(tl); hold(ax,'on');
        for wi=1:2
            p=abs(p_circle_f(:,wi,fi)); p=p/max(p(theta_deg==0));
            C=[.0 .35 .68;.83 .35 .05];
            plot(ax,theta_deg,20*log10(max(p,realmin)),'Color',C(wi,:),'LineWidth',2, ...
                'DisplayName',sprintf('%d cm exact',round(100*cfg.total_widths_m(wi))));
            ii=1:24:numel(theta_deg);
            plot(ax,theta_deg(ii),20*log10(max(abs(D_sinc(ii,wi,fi)),realmin)), ...
                '.','Color',C(wi,:),'MarkerSize',7,'HandleVisibility','off');
        end
        yline(ax,-3,'k:','-3 dB','HandleVisibility','off');
        xline(ax,atan2d(300,100)*[-1 1],'k:','HandleVisibility','off');
        grid(ax,'on'); xlim(ax,[-90 90]); ylim(ax,[-45 1]);
        title(ax,sprintf('%d kHz',cfg.frequencies_hz(fi)/1000));
        xlabel(ax,'Angle from upward (deg)'); ylabel(ax,'Axis-normalized amplitude (dB)');
        if fi==1, legend(ax,'Location','south','FontSize',9); end
    end
    savefig(fig,fullfile(out,'directivity.fig'));
    exportgraphics(fig,fullfile(out,'directivity.png'),'Resolution',180); close(fig);
    receipt=struct('run_id',run_id,'status','FIGURES_AND_DATA_COMPLETE_AUDIT_PENDING', ...
        'elapsed_s',toc(t0),'failed_checks',0,'result_dir',out, ...
        'clim_dB_re_1Pa',clim_pa,'quadrature_error',quad_error, ...
        'farfield_sinc_max_abs_error',farfield_error,'half_3db_deg',half_3db_deg, ...
        'first_null_deg',first_null_deg,'matlab_release',version);
    fid=fopen(fullfile(out,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
    fprintf('%s\n',jsonencode(receipt));
catch ME
    fid=fopen(fullfile(out,'error.txt'),'w'); fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    rethrow(ME);
end
    function format_axes(ax)
        set(ax,'YDir','reverse','FontSize',10,'Layer','top');
        xlim(ax,[-300 300]); ylim(ax,[0 100]);
        xlabel(ax,'Horizontal x (m)'); ylabel(ax,'Depth z (m), positive downward');
    end
    function [xe,ze,clipped]=ray_endpoint(th)
        s_top=100/cosd(th); s_side=300/max(abs(sind(th)),eps);
        s=min(s_top,s_side); xe=s*sind(th); ze=100-s*cosd(th);
        clipped=s_side<s_top;
    end
end
