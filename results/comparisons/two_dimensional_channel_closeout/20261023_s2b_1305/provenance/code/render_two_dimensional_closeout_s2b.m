function render_two_dimensional_closeout_s2b(matpath,out)
% Saved-data-only F01-F07 renderer. No solver, interpolation or field calls.
assert(~isfolder(out),'Refusing overwrite'); mkdir(out);
t=load(matpath,'plot_data'); d=t.plot_data; f=d.line; U=d.U; g=d.grid;
names={'PE-mode 0','BIE','FullFK'}; styles={'-','--',':'};
cols=[0 .35 .65;.85 .32 .10;.1 .55 .3]; labels={'flat','weak A=0.01 m','medium A=0.05 m','strong A=0.20 m'};
fig=figure('Visible','off','Position',[50 50 1200 620],'Color','w'); tiledlayout(fig,1,2);
ax=nexttile; hold(ax,'on'); sx=linspace(-60,60,1201)';
for j=2:4, eta=closeout_s2b_surface(sx,[0 .01 .05 .2]*(1:4==j)'); plot(ax,sx,eta,'DisplayName',labels{j}); end
plot(ax,0,100,'kp','MarkerSize',12,'DisplayName','Tx (0,100 m)'); plot(ax,g.x_m,g.z_m,'ko','DisplayName','Rx 9 true points');
set(ax,'YDir','reverse'); ylim(ax,[-3 105]); xlabel(ax,'x (m)'); ylabel(ax,'depth z (m), positive down'); title(ax,'Frozen 2-D physical geometry'); legend(ax,'Location','southwest'); grid(ax,'on');
ax=nexttile; scatter(ax,g.x_m,g.z_m,50,'k','filled'); text(ax,g.x_m+.008,g.z_m,string(g.receiver_id));
set(ax,'YDir','reverse'); xlim(ax,[-.25 .3]); ylim(ax,[2.83 3.17]); xlabel(ax,'x (m)'); ylabel(ax,'depth z (m)'); title(ax,'P11-P33: no shifted coordinates'); grid(ax,'on'); save_owned(fig,out,'F01_geometry');
fig=figure('Visible','off','Position',[50 50 1200 720],'Color','w'); tiledlayout(fig,2,2);
for m=1:3
    a=nexttile(1); hold(a,'on'); plot(a,f.x_m,abs(f.(method(m)){1}),styles{m},'Color',cols(m,:),'DisplayName',names{m});
    a=nexttile(2); hold(a,'on'); plot(a,f.x_m,angle(f.(method(m)){1}),styles{m},'Color',cols(m,:));
    a=nexttile(3); hold(a,'on'); plot(a,f.x_m,abs(f.(method(m)){1}-f.BIE{1}),styles{m},'Color',cols(m,:));
    a=nexttile(4); hold(a,'on'); plot(a,f.x_m,angle(f.(method(m)){1}.*conj(f.BIE{1})),styles{m},'Color',cols(m,:));
end
titles={'Flat reflected magnitude (unit-source response)','Flat native wrapped phase (rad)','Absolute difference from BIE (same units)','Phase difference from BIE (rad)'};
for j=1:4, ax=nexttile(j); title(ax,titles{j}); xlabel(ax,'x (m), z=3 m; M99 411 real samples'); grid(ax,'on'); end
legend(nexttile(1),'Location','best'); save_owned(fig,out,'F02_flat');
fig=figure('Visible','off','Position',[50 50 1250 700],'Color','w'); tiledlayout(fig,2,2);
lim=max([abs(d.strong.PE(:));abs(d.strong.BIE(:))]);
for m=1:2
    a=nexttile(m); imagesc(a,d.strong.x_m,d.strong.z_m,abs(d.strong.(method(m))).'); set(a,'YDir','reverse'); clim(a,[0 lim]); colorbar(a); title(a,[names{m} ' strong: |u|, unit-source response']); xlabel(a,'x (m)'); ylabel(a,'depth z (m)');
    a=nexttile(m+2); imagesc(a,d.strong.x_m,d.strong.z_m,angle(d.strong.(method(m))).'); set(a,'YDir','reverse'); clim(a,[-pi pi]); colorbar(a); title(a,[names{m} ' strong: native wrapped phase (rad)']); xlabel(a,'x (m)'); ylabel(a,'depth z (m)');
end
sgtitle(fig,'Saved true 411 x 45 field; FullFK continuous field absent (no interpolation)'); save_owned(fig,out,'F03_saved_strong_xz');
fig=figure('Visible','off','Position',[50 50 1300 740],'Color','w'); tiledlayout(fig,2,3);
for c=2:4
    for m=1:3
        a=nexttile(c-1); hold(a,'on'); plot(a,f.x_m,abs(f.(method(m)){c}),styles{m},'Color',cols(m,:),'DisplayName',names{m});
        a=nexttile(c+2); hold(a,'on'); plot(a,f.x_m,angle(f.(method(m)){c}),styles{m},'Color',cols(m,:));
    end
    a=nexttile(c-1); title(a,labels{c}); ylabel(a,'|u reflected| (unit-source response)'); xlabel(a,'x (m), z=3 m'); grid(a,'on');
    a=nexttile(c+2); ylabel(a,'native wrapped phase (rad)'); xlabel(a,'x (m), z=3 m'); grid(a,'on');
end
legend(nexttile(1),'Location','best'); sgtitle(fig,'Unnormalized same-source native fields: 411 real M99 samples'); save_owned(fig,out,'F04_native_lines');
fig=figure('Visible','off','Position',[50 50 1250 760],'Color','w'); tiledlayout(fig,2,3);
pairs=[1 2;3 2;1 3];
for p=1:3
    E=zeros(4,9); PH=E;
    for c=1:4, a=U(:,c,pairs(p,1)); b=U(:,c,pairs(p,2)); E(c,:)=abs(a-b)./abs(b); PH(c,:)=angle(a.*conj(b)); end
    a=nexttile(p); imagesc(a,E); colorbar(a); title(a,[names{pairs(p,1)} ' - ' names{pairs(p,2)}]); set(a,'XTick',1:9,'XTickLabel',g.receiver_id,'YTick',1:4,'YTickLabel',labels); xtickangle(a,45);
    a=nexttile(p+3); imagesc(a,PH); colorbar(a); title(a,'Signed phase difference (rad)'); set(a,'XTick',1:9,'XTickLabel',g.receiver_id,'YTick',1:4,'YTickLabel',labels); xtickangle(a,45);
end
sgtitle(fig,'Top: relative complex error; bottom: phase (rad). 108 true values; separate color scales'); save_owned(fig,out,'F05_point_residuals');
fig=figure('Visible','off','Position',[50 50 1300 720],'Color','w'); tiledlayout(fig,2,4); five=[2 4 5 6 8];
for c=1:4
    for m=1:3
        a=nexttile(c); hold(a,'on'); plot(a,1:5,abs(U(five,c,m)),[styles{m} 'o'],'Color',cols(m,:),'DisplayName',names{m});
        a=nexttile(c+4); hold(a,'on'); plot(a,1:5,angle(U(five,c,m)),[styles{m} 'o'],'Color',cols(m,:));
    end
    a=nexttile(c); title(a,labels{c}); ylabel(a,'|u| (unit-source response)'); set(a,'XTick',1:5,'XTickLabel',g.receiver_id(five)); grid(a,'on');
    a=nexttile(c+4); ylabel(a,'native phase (rad)'); set(a,'XTick',1:5,'XTickLabel',g.receiver_id(five)); grid(a,'on');
end
legend(nexttile(1),'Location','best'); sgtitle(fig,'Five fixed representative positions; connecting segments are guides only'); save_owned(fig,out,'F06_five_receivers');
fig=figure('Visible','off','Position',[50 50 1300 500],'Color','w'); tiledlayout(fig,1,3);
for p=1:3
    a=nexttile; hold(a,'on');
    line=d.line_metrics; ar=d.array_metrics;
    ca=string(method(pairs(p,1))); rb=string(method(pairs(p,2)));
    s=line.candidate==ca&line.reference==rb; q=line(s&line.weight_kind=="spatial_uniform_411",:); w=line(s&line.weight_kind=="legacy_incident_energy_M99",:);
    y=ar.complex_l2(ar.candidate==ca&ar.reference==rb);
    plot(a,1:4,q.complex_l2,'-o','DisplayName','411 spatial equal weights'); plot(a,1:4,w.complex_l2,'--s','DisplayName','411 legacy incident weights'); plot(a,1:4,y,':^','DisplayName','9-point equal weights');
    set(a,'YScale','log','XTick',1:4,'XTickLabel',{'flat','weak','medium','strong'}); ylabel(a,'relative complex L2'); title(a,[names{pairs(p,1)} ' - ' names{pairs(p,2)}]); grid(a,'on'); legend(a,'Location','best');
end
sgtitle(fig,'Different fixed supports/weights; REFERENCE LIMITED, no physical PASS'); save_owned(fig,out,'F07_metrics');
end

function s=method(m)
names={'PE','BIE','FullFK'}; s=names{m};
end
function save_owned(fig,out,stem)
drawnow; savefig(fig,fullfile(out,[stem '.fig'])); exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',160); close(fig);
end
