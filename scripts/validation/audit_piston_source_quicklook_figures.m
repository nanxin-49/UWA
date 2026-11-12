function receipt = audit_piston_source_quicklook_figures(run_dir)
%AUDIT_PISTON_SOURCE_QUICKLOOK_FIGURES Independently reconstruct FIG encodings.
arguments
    run_dir (1,:) char
end
t0=tic; out=fullfile(run_dir,'figure_audit_delivered'); assert(~isfolder(out)); mkdir(out);
figdir=fullfile(run_dir,'figures_delivered'); checks=[]; names={};
for fk=[4 6 8]
    s=load(fullfile(run_dir,sprintf('field_%dkHz.mat',fk)));
    for other=[3 2]
        stem=sprintf('field_%dkHz_20vs%02dcm',fk,round(100*s.cfg.total_widths_m(other)));
        fig=openfig(fullfile(figdir,[stem '.fig']),'new','invisible');
        ims=findobj(fig,'Type','image');
        expected={20*log10(abs(s.p_field_f(:,:,1))),20*log10(abs(s.p_field_f(:,:,other))), ...
            20*log10(abs(s.p_field_f(:,:,1))./abs(s.p_field_f(:,:,other))), ...
            20*log10(abs(s.p_near_f(:,:,1))),20*log10(abs(s.p_near_f(:,:,other))), ...
            20*log10(abs(s.p_near_f(:,:,1))./abs(s.p_near_f(:,:,other)))};
        assert(numel(ims)==6); matches=false(1,6); err=0; scale_ok=true;
        for ii=1:6
            dd=ims(ii).CData; best=Inf; bestj=0;
            for jj=1:6
                if ~matches(jj) && isequal(size(dd),size(expected{jj}))
                    delta=max(abs(dd(:)-expected{jj}(:)));
                    if delta<best, best=delta; bestj=jj; end
                end
            end
            assert(bestj>0); matches(bestj)=true; err=max(err,best);
            ax=ancestor(ims(ii),'axes');
            if ismember(bestj,[3 6]), limits=[-24 24]; else, limits=[0 70]; end
            scale_ok=scale_ok && isequal(ax.CLim,limits) && strcmp(ax.YDir,'reverse');
        end
        info=imfinfo(fullfile(figdir,[stem '.png']));
        checks=[checks;err scale_ok info.Width info.Height]; names{end+1}=stem; %#ok<AGROW>
        if fk==8 && other==3
            drawnow;
            exportgraphics(fig,fullfile(out,'field_8kHz_reopened.png'),'Resolution',150);
            savefig(fig,fullfile(out,'field_8kHz_reopened.fig'));
        end
        close(fig);
    end
end
% Every data curve on the line comparison must match a fresh MAT transform.
fig=openfig(fullfile(figdir,'receiver_line_comparison.fig'),'new','invisible');
lines=findobj(fig,'Type','line'); expected_count=0; missing_count=0; negative_gain_visible=true;
for fk=[4 6 8]
    s=load(fullfile(run_dir,sprintf('field_%dkHz.mat',fk)));
    for row=1:3
        for wi=1:3
            if row==3 && wi==1, continue; end
            if row==1
                yy=abs(s.p_rx_f(:,wi));
            elseif row==2
                yy=20*log10(abs(s.p_rx_f(:,wi))/abs(s.p_rx_f(s.x_rx_m==0,wi)));
            else
                yy=20*log10(abs(s.p_rx_f(:,1))./abs(s.p_rx_f(:,wi)));
            end
            matched=false;
            for ii=1:numel(lines)
                xx=lines(ii).XData; y=lines(ii).YData;
                if isequal(xx(:),s.x_rx_m) && isequal(y(:),yy)
                    matched=true;
                    if fk==8 && row==3
                        ax=ancestor(lines(ii),'axes'); negative_gain_visible=negative_gain_visible && ax.YLim(1)<=min(yy);
                    end
                    break;
                end
            end
            expected_count=expected_count+1; missing_count=missing_count+~matched;
        end
    end
end
close(fig);
receipt=struct('status','FIG_DATA_RECONSTRUCTION_COMPLETE','elapsed_s',toc(t0), ...
    'maps_checked',{names},'max_CData_error_dB',max(checks(:,1)), ...
    'common_color_scales_and_depth_ok',all(checks(:,2)), ...
    'receiver_curves_checked',expected_count,'receiver_curves_missing',missing_count, ...
    'negative_offaxis_gain_visible',negative_gain_visible, ...
    'failed_checks',nnz(checks(:,1)>1e-12)+nnz(~checks(:,2))+missing_count+~negative_gain_visible, ...
    'PNG_visual_review','separate human-visible model review recorded in visual_review.md');
fid=fopen(fullfile(out,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt)); assert(receipt.failed_checks==0);
end
