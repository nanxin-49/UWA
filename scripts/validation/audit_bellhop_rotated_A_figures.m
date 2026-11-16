function receipt=audit_bellhop_rotated_A_figures(run_dir,figure_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_A_FIGURES Bind complete native vectors and raw colors.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic;
s=load(fullfile(run_dir,'assembly','result.mat'),'total','direct','surface_family','bottom_family','envelope','plan');
a=load(fullfile(run_dir,'assembly','sum_through_order12.mat'),'total'); x=s.plan.x_m(:); z=s.plan.depths_m(:);
checks=table(strings(0,1),false(0,1),zeros(0,1),'VariableNames',{'check','pass','max_difference'});
files=dir(fullfile(figure_dir,'F*.fig'));
for j=1:numel(files)
    path=fullfile(files(j).folder,files(j).name); fig=openfig(path,'new','invisible'); clean=onCleanup(@()close(fig));
    images=findall(fig,'Type','image');
    for k=1:numel(images)
        u=images(k).UserData; rows=u.rows; field=u.field;
        if strcmp(field,'tail_near'), expected=abs(s.total(rows,:)-a.total(rows,:))./s.envelope(rows,:);
        else, if strcmp(field,'total_near'), field='total'; end; expected=-20*log10(abs(s.(field)(rows,:))); end
        name=[files(j).name '_' images(k).Tag];
        add([name '_raw_colors'],isequaln(images(k).CData,expected),maximum(images(k).CData,expected));
        add([name '_all_original_x'],isequal(images(k).XData(:),x),max(abs(images(k).XData(:)-x)));
        add([name '_all_original_z'],isequal(images(k).YData(:),z(rows)),max(abs(images(k).YData(:)-z(rows))));
        add([name '_source_and_missing_alpha'],isequal(logical(images(k).AlphaData),isfinite(expected)),0);
        add([name '_uniform_strip_no_hidden_resample'],max(abs(diff(z(rows))-mean(diff(z(rows)))))<1e-9,0);
    end
    if isempty(images)
        lines=findall(fig,'Type','line'); ix=find(x==0);
        for k=1:numel(lines)
            name=lines(k).DisplayName;
            if strcmp(name,'N12'), expected=abs(a.total(:,ix));
            elseif strcmp(name,'N16'), expected=abs(s.total(:,ix));
            else, expected=angle(a.total(:,ix).*conj(s.total(:,ix))); expected(abs(s.total(:,ix))<.05*s.envelope(:,ix))=NaN;
            end
            add(sprintf('%s_line%d_native_z',files(j).name,k),isequal(lines(k).XData(:),z),0);
            add(sprintf('%s_line%d_raw_y',files(j).name,k),isequaln(lines(k).YData(:),expected),maximum(lines(k).YData(:),expected));
        end
    end
    [~,stem]=fileparts(path); png=fullfile(figure_dir,[stem '.png']); info=imfinfo(png); add([stem '_PNG_shape'],info.Width>1000&&info.Height>500,0);
    clear clean
end
receipt=struct('status','SAVED_FIGURE_FULL_VECTOR_AUDIT','failed_checks',nnz(~checks.pass),'checks_count',height(checks), ...
    'elapsed_s',toc(start),'figure_dir',figure_dir,'output_dir',output_dir,'source_MAT_SHA256',black_file_hash(fullfile(run_dir,'assembly','result.mat')), ...
    'preserved_failure','Original render checkpoint has12 false failures because endpoints were compared with stored full vectors; no numerical data changed', ...
    'physics_solve_count',0); writetable(checks,fullfile(output_dir,'checks.csv'));
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,pass,diff)
        checks(end+1,:)={string(name),logical(pass),diff};
    end
end
function v=maximum(a,b)
use=isfinite(a)&isfinite(b); v=max(abs(a(use)-b(use))); if isempty(v), v=0; end
end