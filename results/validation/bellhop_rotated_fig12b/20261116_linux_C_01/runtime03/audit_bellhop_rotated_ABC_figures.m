function receipt=audit_bellhop_rotated_ABC_figures(A_file,B_file,C_file,figure_dir,output_dir)
% Independent raw MAT to reopened FIG audit; no propagation.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic; F={load(A_file,'total','plan'),load(B_file,'total','plan'),load(C_file,'total','envelope','plan')};
x=F{3}.plan.x_m(:); z=F{3}.plan.depths_m(:); checks=table(strings(0,1),false(0,1),zeros(0,1),'VariableNames',{'check','pass','max_difference'});
files=dir(fullfile(figure_dir,'F*.fig')); assert(numel(files)==3);
for j=1:3
 fig=openfig(fullfile(files(j).folder,files(j).name),'new','invisible'); cleanup=onCleanup(@()close(fig)); ims=findall(fig,'Type','image');
 for k=1:numel(ims)
  u=ims(k).UserData; r=u.rows; expected=-20*log10(abs(F{u.group}.total(r,:))); name=sprintf('%s_image%d',files(j).name,k);
  add([name '_colors'],isequaln(ims(k).CData,expected),difference(ims(k).CData,expected)); add([name '_all_x'],isequal(ims(k).XData(:),x),0); add([name '_all_z'],isequal(ims(k).YData(:),z(r)),0);
  add([name '_alpha'],isequal(logical(ims(k).AlphaData),isfinite(expected)),0); add([name '_CLim'],isequal(ims(k).Parent.CLim,[20 70]),0);
 end
 if isempty(ims)
  ls=findall(fig,'Type','line'); ix=find(x==0); assert(numel(ls)==4);
  for k=1:numel(ls)
   u=ls(k).UserData;
   if strcmp(u.kind,'amplitude'), expected=abs(F{u.group}.total(:,ix)); else, expected=angle(F{2}.total(:,ix).*conj(F{3}.total(:,ix))); expected(abs(F{3}.total(:,ix))<.05*F{3}.envelope(:,ix))=NaN; end
   add(sprintf('line%d_all_z',k),isequal(ls(k).XData(:),z),0); add(sprintf('line%d_values',k),isequaln(ls(k).YData(:),expected),difference(ls(k).YData(:),expected));
  end
 end
 [~,stem]=fileparts(files(j).name); info=imfinfo(fullfile(figure_dir,[stem '.png'])); add([stem '_PNG'],info.Width>1000&&info.Height>500,0); clear cleanup
end
receipt=struct('status','ABC_FIGURE_INDEPENDENT_AUDIT','failed_checks',nnz(~checks.pass),'checks_count',height(checks),'elapsed_s',toc(start),'figure_dir',figure_dir,'new_native_processes',0);
writetable(checks,fullfile(output_dir,'checks.csv')); fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
 function add(name,pass,v), checks(end+1,:)={string(name),logical(pass),v}; end
end
function v=difference(a,b)
use=isfinite(a)&isfinite(b); v=max(abs(a(use)-b(use))); if isempty(v), v=0; end
end
