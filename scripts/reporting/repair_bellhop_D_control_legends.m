function receipt=repair_bellhop_D_control_legends(project_root,out)
%REPAIR_BELLHOP_D_CONTROL_LEGENDS Saved data only; preserve every original.
t=tic; assert(~isfolder(out));mkdir(out);addpath(fullfile(project_root,'scripts','validation','support'));
d201=fullfile(project_root,'results','validation','bellhop_1500m_applicability','20261118_D4000_C_retry02');d801=fullfile(project_root,'results','validation','bellhop_1500m_applicability','20261119_D4000_C801_01','returned');
a=load(fullfile(d201,'independent_audit.mat'),'fields','ray_end'); b=load(fullfile(d801,'comparison.mat'),'p201','p401','p801','ref','E','x','metrics','gridtable');
old={fullfile(d201,'D4000_C_diagnostic'),fullfile(d801,'F01_complex_error_201_401_801'),fullfile(d801,'F02_angular_error_scaling')}; source=table(strings(6,1),strings(6,1),'VariableNames',{'path','sha256'});n=0;
for k=1:3,for ext={'.fig','.png'},n=n+1;source.path(n)=[old{k} ext{1}];source.sha256(n)=black_file_hash(char(source.path(n)));end;end
expected={}; saved={}; legends=0; curves=0;
f=ownedfig('Saved D201 direct control: legend repair',1200,840); ax=panels(f); x=a.fields.x_m; ref=a.fields.reference_p_real+1i*a.fields.reference_p_imag;p=a.fields.Bellhop_p_real+1i*a.fields.Bellhop_p_imag;
h=plot(ax(1),x,abs(ref),'k-o',x,abs(p),'b--s');set(h(1),'Tag','analytic_D');set(h(2),'Tag','native_D201');ylabel(ax(1),'Pressure / 1 m reference');
lg=legend(ax(1),h,{'Analytic D','Bellhop 201'},'Orientation','horizontal','AutoUpdate','off');manuallegend(lg,[.28 .012 .44 .035]);legends=legends+1;
h=plot(ax(2),x,a.fields.complex_difference_over_E,'b-o');set(h,'Tag','D201_complex_E');yline(ax(2),.002,'r:','Frozen max gate','HandleVisibility','off');ylabel(ax(2),'|complex difference| / E');
h=plot(ax(3),x,a.fields.phase_difference_rad,'b-o');set(h,'Tag','D201_phase');yline(ax(3),.01,'r:','HandleVisibility','off');yline(ax(3),-.01,'r:','HandleVisibility','off');ylabel(ax(3),'Phase difference (rad)');
h=plot(ax(4),a.ray_end.launch_angle_deg,a.ray_end.end_tau_error_s*1e6,'k-');set(h,'Tag','D201_time');xlabel(ax(4),'Native launch angle (deg)');ylabel(ax(4),'Terminal time error (microseconds)');
for j=1:3,xlabel(ax(j),'Independent link offset x (m)');xlim(ax(j),[-5 205]);end
exp={{'analytic_D',x,abs(ref)},{'native_D201',x,abs(p)},{'D201_complex_E',x,a.fields.complex_difference_over_E},{'D201_phase',x,a.fields.phase_difference_rad},{'D201_time',a.ray_end.launch_angle_deg,a.ray_end.end_tau_error_s*1e6}};
curves=curves+saveaudit(f,out,'R01_D201_legend_fixed',exp);saved{end+1}='R01_D201_legend_fixed';expected{end+1}=exp;
P=[b.p201(:),b.p401(:),b.p801(:)]; ref=b.ref(:); E=b.E(:);x=b.x(:);beams=[201 401 801];
f=ownedfig('U0 direct 4000 Hz: saved 201 / 401 / 801',1200,840);ax=panels(f);vals={real(P-ref)./E,imag(P-ref)./E,abs(P-ref)./E,angle(P.*conj(ref))};labels={'Re(complex residual) / E','Im(complex residual) / E','|complex residual| / E','Phase error (rad)'};exp={};
for j=1:4
 hh=plot(ax(j),x,vals{j}(:,1),'k--o',x,vals{j}(:,2),'b-.s',x,vals{j}(:,3),'r-^');for k=1:3,tag=sprintf('B%d_P%d',beams(k),j);set(hh(k),'Tag',tag);exp{end+1}={tag,x,vals{j}(:,k)};end
 xlabel(ax(j),'Independent link offset x (m)');ylabel(ax(j),labels{j});xlim(ax(j),[-5 205]);if j==1,lg=legend(ax(j),hh,{'201 saved','401 saved','801 saved'},'Orientation','horizontal','AutoUpdate','off');manuallegend(lg,[.26 .012 .48 .035]);legends=legends+1;end
 if j==3,yline(ax(j),.002,'k:','Frozen max gate','HandleVisibility','off');end
end
curves=curves+saveaudit(f,out,'R02_complex_error_legend_fixed',exp);saved{end+1}='R02_complex_error_legend_fixed';expected{end+1}=exp;
f=ownedfig('Frozen angular grids: error scale and phase prediction',1280,560);ax(1)=axes('Parent',f,'Position',[.075 .22 .39 .68]);ax(2)=axes('Parent',f,'Position',[.575 .22 .39 .68]);
hh=loglog(ax(1),b.metrics.delta_angle_rad,b.metrics.complex_L2,'k-o',b.metrics.delta_angle_rad,b.metrics.max_complex_over_E,'b-s',b.metrics.delta_angle_rad,b.metrics.complex_L2(1)*(b.metrics.delta_angle_rad/b.metrics.delta_angle_rad(1)).^2,'r--');
ys={b.metrics.complex_L2,b.metrics.max_complex_over_E,b.metrics.complex_L2(1)*(b.metrics.delta_angle_rad/b.metrics.delta_angle_rad(1)).^2};tags={'complex_L2','max_E','h2_guide'};exp={};for k=1:3,set(hh(k),'Tag',tags{k});exp{end+1}={tags{k},b.metrics.delta_angle_rad,ys{k}};end
xlabel(ax(1),'Angular interval h (rad)');ylabel(ax(1),'Dimensionless complex error');lg=legend(ax(1),hh,{'Complex L2','Max / E','h^2 guide (201 anchor)'},'Orientation','horizontal','AutoUpdate','off');manuallegend(lg,[.045 .025 .48 .055]);legends=legends+1;
pred=reshape(b.gridtable.second_order_phase_prediction_rad,5,3);ph=reshape(b.gridtable.measured_phase_error_rad,5,3);hh=plot(ax(2),pred(:,1),ph(:,1),'ko',pred(:,2),ph(:,2),'bs',pred(:,3),ph(:,3),'r^');hold(ax(2),'on');lims=[min(pred,[],'all'),0];eq=plot(ax(2),lims,lims,'k:');
for k=1:3,tag=sprintf('prediction_%d',beams(k));set(hh(k),'Tag',tag);exp{end+1}={tag,pred(:,k),ph(:,k)};end;set(eq,'Tag','equality');exp{end+1}={'equality',lims,lims};
xlabel(ax(2),'Second-order phase prediction (rad)');ylabel(ax(2),'Measured phase error (rad)');lg=legend(ax(2),[hh;eq],{'201','401','801','Equality'},'Orientation','horizontal','AutoUpdate','off');manuallegend(lg,[.575 .025 .39 .055]);legends=legends+1;
for j=1:2,grid(ax(j),'on');set(ax(j),'FontSize',10);end
curves=curves+saveaudit(f,out,'R03_angular_scaling_legend_fixed',exp);saved{end+1}='R03_angular_scaling_legend_fixed';expected{end+1}=exp;
for k=1:height(source),assert(strcmpi(black_file_hash(char(source.path(k))),source.sha256(k)),'Original changed');end
writetable(source,fullfile(out,'original_figure_sha256.csv'));receipt=struct('status','SAVED_FIGURE_REPAIR_DATA_AUDIT_PASS','new_acoustic_solves',0,'original_files_checked',height(source),'original_files_changed',0,'repaired_figures',3,'exact_curve_checks',curves,'legend_checks',legends,'all_legends_inside_canvas',true,'legend_AutoUpdate','off','visual_review','PENDING','elapsed_s',toc(t),'matlab_version',version);
save(fullfile(out,'repair_audit.mat'),'receipt','source','saved','expected');fid=fopen(fullfile(out,'repair_receipt.json'),'w');fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));fclose(fid);
fprintf('FIGURE_REPAIR figures3 curves%d legends%d originalSHA6PASS solves0 elapsed%.3fs\n',curves,legends,toc(t));
end
function f=ownedfig(titletext,w,h)
f=figure('Visible','off','Color','w','Position',[50 50 w h]);annotation(f,'textbox',[0 .95 1 .04],'String',titletext,'HorizontalAlignment','center','EdgeColor','none','FontSize',14,'Interpreter','none');
end
function ax=panels(f)
pos={[.08 .57 .39 .34],[.57 .57 .39 .34],[.08 .13 .39 .34],[.57 .13 .39 .34]};ax=gobjects(1,4);for k=1:4,ax(k)=axes('Parent',f,'Position',pos{k},'FontSize',10);grid(ax(k),'on');end
end
function manuallegend(lg,pos)
set(lg,'Units','normalized','Position',pos,'AutoUpdate','off','FontSize',9);assert(pos(1)>=0&&pos(2)>=0&&pos(1)+pos(3)<=1&&pos(2)+pos(4)<=1);
end
function n=saveaudit(f,out,stem,exp)
drawnow;savefig(f,fullfile(out,[stem '.fig']));exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',140);close(f);
f=openfig(fullfile(out,[stem '.fig']),'invisible');for k=1:numel(exp),e=exp{k};o=findobj(f,'Tag',e{1});assert(isscalar(o)&&isequal(o.XData,e{2}(:)')&&isequal(o.YData,e{3}(:)'),'Curve mismatch');end
lg=findall(f,'Type','legend');for k=1:numel(lg),assert(strcmp(lg(k).AutoUpdate,'off'));p=lg(k).Position;assert(all(p(1:2)>=0)&&p(1)+p(3)<=1&&p(2)+p(4)<=1);assert(~any(startsWith(string(lg(k).String),'data')));end;close(f);n=numel(exp);
end