function out = render_two_dimensional_closeout_s3a_three_seed(root, dest)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S3A_THREE_SEED Saved-field comparison only.
arguments
    root (1,:) char
    dest (1,:) char
end
assert(isfolder(root) && ~isfolder(dest));
parent=fullfile(root,'results','comparisons','two_dimensional_channel_closeout');
runs={'20261025_s3a_seed260001_01','20261025_s3a_seed260002_01', ...
    '20261026_s3a_seed260003_01'};
seeds=[260001 260002 260003];
support=["nine_9" "line_411"];
rows=table;
geometry=table;
for k=1:3
    run=fullfile(parent,runs{k});
    a=jsondecode(fileread(fullfile(run,'audit','independent_saved_data_audit_v3.json')));
    assert(strcmp(a.status,sprintf('S3A_SEED%d_SAVED_DATA_AUDIT_PASS',seeds(k))));
    fig_audit=fullfile(run,'audit','independent_figure_content.json');
    if isfile(fig_audit)
        f=jsondecode(fileread(fig_audit));
        assert(strcmp(f.status,sprintf('S3A_SEED%d_FIGURE_CONTENT_PASS',seeds(k))));
    end
    assert(all(a.fig_reexport_pixel_differences==0));
    m=readtable(fullfile(run,'analysis_v2','model_metrics.csv'),TextType='string');
    b=readtable(fullfile(run,'analysis_v2','BIE_physical_response.csv'),TextType='string');
    assert(height(m)==4 && height(b)==2);
    for j=1:2
        sid=support(j);
        pe=m(m.method=="PE" & m.support_id==sid,:);
        fk=m(m.method=="FullFK" & m.support_id==sid,:);
        pb=b(b.support_id==sid,:);
        assert(height(pe)==1 && height(fk)==1 && height(pb)==1);
        assert(pe.above_5x_empirical_resolution && fk.above_5x_empirical_resolution ...
            && pb.above_5x_empirical_resolution);
        rows=[rows;table(seeds(k),sid,pb.n,pb.BIE_physical_response_L2, ...
            pb.BIE_level_resolution_relative,pe.complex_L2,pe.amplitude_L2, ...
            pe.phase_L2,pe.phase_RMS_rad,pe.total_empirical_resolution, ...
            fk.complex_L2,fk.amplitude_L2,fk.phase_L2,fk.phase_RMS_rad, ...
            fk.total_empirical_resolution, ...
            'VariableNames',{'seed','support_id','n','BIE_physical_response_L2', ...
            'BIE_reference_resolution','PE_complex_L2','PE_amplitude_L2', ...
            'PE_phase_L2','PE_phase_RMS_rad','PE_numeric_resolution', ...
            'FullFK_complex_L2','FullFK_amplitude_L2','FullFK_phase_L2', ...
            'FullFK_phase_RMS_rad','FullFK_numeric_resolution'})]; %#ok<AGROW>
    end
    p=readtable(fullfile(run,'input',sprintf('seed_%d_coefficients.csv',seeds(k))));
    xx=(-80:.001:80)';
    [eta,sl,cv,raw]=closeout_s3a_pm_surface(p,xx);
    geometry=[geometry;table(seeds(k),max(abs(raw)),max(abs(eta)), ...
        max(abs(sl)),max(abs(cv)),min(eta), ...
        'VariableNames',{'seed','raw_max_abs_height_m', ...
        'embedded_max_abs_height_m','embedded_max_abs_slope', ...
        'embedded_max_abs_curvature_per_m','embedded_min_height_m'})]; %#ok<AGROW>
end
mkdir(dest);
writetable(rows,fullfile(dest,'three_seed_comparison.csv'));
writetable(geometry,fullfile(dest,'three_seed_geometry.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1250 520]);
t=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
names={'BIE PM / 同源平面','PE / BIE','完整 FullFK / BIE'};
cols={'BIE_physical_response_L2','PE_complex_L2','FullFK_complex_L2'};
for q=1:3
    ax=nexttile(t); hold(ax,'on');
    for j=1:2
        v=zeros(3,1);
        for k=1:3
            z=rows(rows.seed==seeds(k) & rows.support_id==support(j),:);
            v(k)=z.(cols{q});
        end
        if j==1
            plot(ax,1:3,v,'-o','LineWidth',1.6,'MarkerSize',8, ...
                'DisplayName','真实九点');
        else
            plot(ax,1:3,v,'-s','LineWidth',1.6,'MarkerSize',8, ...
                'DisplayName','411点线');
        end
    end
    if q>1, set(ax,'YScale','log'); end
    xlim(ax,[.8 3.2]); xticks(ax,1:3); xticklabels(ax,{'260001','260002','260003'});
    xlabel(ax,'固定Fourier系数 seed'); ylabel(ax,'独立分母复 L2');
    title(ax,names{q}); grid(ax,'on');
    if q==1, legend(ax,'Location','best'); end
end
title(t,'S3A：三个固定PM实现，同一海面两种空间支持');
figfile=fullfile(dest,'F08_S3A_three_seed_summary.fig');
pngfile=fullfile(dest,'F08_S3A_three_seed_summary.png');
savefig(fig,figfile); exportgraphics(fig,pngfile,'Resolution',160); close(fig);
out=struct('status','THREE_SEED_SAVED_DATA_SUMMARY_PENDING_INDEPENDENT_AUDIT', ...
    'rows',height(rows),'seed_count',3,'figure',pngfile);
end
