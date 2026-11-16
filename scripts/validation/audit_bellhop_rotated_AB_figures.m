function receipt=audit_bellhop_rotated_AB_figures(A_file,B_file,figure_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_AB_FIGURES Complete native vectors, colors and curves.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic;
A=load(A_file,'total','plan'); B=load(B_file,'total','envelope','plan'); x=B.plan.x_m(:); z=B.plan.depths_m(:);
assert(isequal(A.plan.depths_m,B.plan.depths_m)&&isequal(A.plan.x_m,B.plan.x_m));
checks=table(strings(0,1),false(0,1),zeros(0,1),'VariableNames',{'check','pass','max_difference'});
files=dir(fullfile(figure_dir,'F*.fig')); assert(numel(files)==3);
for j=1:numel(files)
    fig=openfig(fullfile(files(j).folder,files(j).name),'new','invisible'); cleanup=onCleanup(@()close(fig)); images=findall(fig,'Type','image');
    for k=1:numel(images)
        u=images(k).UserData; r=u.rows;
        if strcmp(u.field,'A'), expected=-20*log10(abs(A.total(r,:))); limits=[20 70];
        elseif strcmp(u.field,'B'), expected=-20*log10(abs(B.total(r,:))); limits=[20 70];
        else, expected=abs(A.total(r,:)-B.total(r,:))./B.envelope(r,:); limits=[0 .02]; end
        name=[files(j).name '_' images(k).Tag];
        add([name '_raw_colors'],isequaln(images(k).CData,expected),maximum(images(k).CData,expected));
        add([name '_all_x'],isequal(images(k).XData(:),x),0); add([name '_all_z'],isequal(images(k).YData(:),z(r)),0);
        add([name '_missing_alpha'],isequal(logical(images(k).AlphaData),isfinite(expected)),0);
        add([name '_color_limits'],isequal(images(k).Parent.CLim,limits),0);
    end
    if isempty(images)
        lines=findall(fig,'Type','line'); ix=find(x==0);
        for k=1:numel(lines)
            label=lines(k).DisplayName;
            if strcmp(label,'A Q1.25'), expected=abs(A.total(:,ix)); elseif strcmp(label,'B Q0.625'), expected=abs(B.total(:,ix));
            else, expected=angle(A.total(:,ix).*conj(B.total(:,ix))); expected(abs(B.total(:,ix))<.05*B.envelope(:,ix))=NaN; end
            add(sprintf('%s_line%d_z',files(j).name,k),isequal(lines(k).XData(:),z),0);
            add(sprintf('%s_line%d_y',files(j).name,k),isequaln(lines(k).YData(:),expected),maximum(lines(k).YData(:),expected));
        end
    end
    [~,stem]=fileparts(files(j).name); info=imfinfo(fullfile(figure_dir,[stem '.png'])); add([stem '_PNG'],info.Width>1000&&info.Height>500,0); clear cleanup
end
receipt=struct('status','AB_FIGURE_INDEPENDENT_AUDIT','failed_checks',nnz(~checks.pass),'checks_count',height(checks),'elapsed_s',toc(start), ...
    'A_sha256',black_file_hash(A_file),'B_sha256',black_file_hash(B_file),'figure_dir',figure_dir,'output_dir',output_dir,'new_native_processes',0);
writetable(checks,fullfile(output_dir,'checks.csv')); fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,pass,difference)
        checks(end+1,:)={string(name),logical(pass),difference};
    end
end
function v=maximum(a,b)
use=isfinite(a)&isfinite(b); v=max(abs(a(use)-b(use))); if isempty(v), v=0; end
end