function receipt = audit_bellhop_rotated_preflight_figures(figure_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_PREFLIGHT_FIGURES Reopen FIG and compare saved raw data.
arguments
    figure_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
diary(fullfile(output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
data=load(fullfile(figure_dir,'source_data.mat'),'ray_data','axis_data');
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
stems={'F01_rotation_direct_rays','F02_axis_travel_time'};
for k=1:2
    fig=openfig(fullfile(figure_dir,[stems{k} '.fig']),'new','invisible');
    close_owned=onCleanup(@()close(fig));
    axes=findall(fig,'Type','axes');
    if k==1
        add('rotation_axes_count',abs(numel(axes)-3),0);
        for a=1:numel(axes)
            label=axes(a).XLabel.String;
            physical=contains(label,'physical x');
            downward=contains(label,'downward chart');
            id=1+downward;
            run=data.ray_data.runs{id};
            for j=1:3
                rq=run.rays(j).rq; limit=run.cfg.range_limit_m; use=rq(:,1)<=limit;
                r=rq(use,1); q=rq(use,2)-run.cfg.transverse_origin_m;
                name=sprintf('%.2g deg',run.rays(j).angle_deg);
                xx=r; yy=q;
                if physical, xx=q; yy=4950-r; name=['up ' name]; end
                line=findobj(axes(a),'Type','line','DisplayName',name);
                assert(isscalar(line),'Missing/duplicate ray graphics.');
                add(sprintf('F01_axes%d_ray%d_x',a,j),max(abs(line.XData(:)-xx(:))),0);
                add(sprintf('F01_axes%d_ray%d_y',a,j),max(abs(line.YData(:)-yy(:))),0);
            end
        end
    else
        add('time_axes_count',abs(numel(axes)-2),0);
        for a=1:numel(axes)
            id=1+contains(axes(a).Title.String,'Downward');
            for j=1:3
                n=(id-1)*3+j; cfg=data.axis_data.configs{n}; row=data.axis_data.runs{n};
                name=sprintf('Q %.1f m / step %.1f m',cfg.ssp_spacing_m,cfg.step_m);
                line=findobj(axes(a),'Type','line','DisplayName',name);
                assert(isscalar(line));
                expected=abs(row.axis_delay_s-row.reference_delay_s)*1e6;
                add(sprintf('F02_axes%d_curve%d_x',a,j),max(abs(line.XData(:)-cfg.receiver_range_m(:))),0);
                add(sprintf('F02_axes%d_curve%d_y',a,j),max(abs(line.YData(:)-expected(:))),0);
            end
        end
    end
    reconstructed=fullfile(output_dir,[stems{k} '_from_fig.png']);
    exportgraphics(fig,reconstructed,'Resolution',120);
    original=imread(fullfile(figure_dir,[stems{k} '.png'])); replay=imread(reconstructed);
    assert(isequal(size(original),size(replay)));
    add([stems{k} '_pixel_max_difference'],max(abs(double(original(:))-double(replay(:)))),0);
    clear close_owned
end
add('no_unrequested_PDF',numel(dir(fullfile(figure_dir,'*.pdf'))),0);
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status','COMPLETE','failed_checks',nnz(~checks.pass),'elapsed_s',toc(t), ...
    'source_data_sha256',black_file_hash(fullfile(figure_dir,'source_data.mat')), ...
    'scope','Reopened FIG data vs saved raw rays/times, and PNG re-export; not acoustic field qualification', ...
    'result_dir',output_dir,'matlab_release',version);
save(fullfile(output_dir,'result.mat'),'checks','receipt');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end