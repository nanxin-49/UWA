function output=diagnose_deepsea_bie_windows()
%DIAGNOSE_DEEPSEA_BIE_WINDOWS Fixed-resolution windows and actual graph trace.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
outdir=fullfile(root,'results','validation','deepsea_bie_reconstruction_20261014');
file=fullfile(outdir,'bie_window_diagnosis.mat');assert(~isfile(file),'BIE:Exists','No overwrite.');
timer=tic;k=2*pi*6000/1500;sigma=.7;targets=[0 0 .3;.5 .2 1;0 0 5];
truth=deepsea_bie_gaussian_reference(targets,k,-1,sigma,512);windows=struct([]);
for L=[2 3 4]
    h=.25;[z,w]=deepsea_bie_gauss(12);centers=-L+h/2:h:L-h/2;
    x=reshape(centers(:)+z.'*h/2,[],1);wx=repmat(w.'*h/2,numel(centers),1);wx=wx(:);
    [xx,yy]=meshgrid(x,x);ps=[xx(:),yy(:),zeros(numel(xx),1)];
    density=exp(-(ps(:,1).^2+ps(:,2).^2)/(2*sigma^2));
    [S,D]=deepsea_bie_layers(targets,ps,repmat([0 0 1],size(ps,1),1), ...
        reshape(wx*wx.',[],1),density,k,-1);value=D-1i*k*S;
    row=struct('L_m',L,'h_m',h,'quadrature_degree',12,'source_count',numel(density), ...
        'known_density_field_error',norm(value-truth)/norm(truth), ...
        'flat_exterior_correction_norm',norm(truth-value)/norm(truth));
    if isempty(windows),windows=row;else,windows(end+1)=row;end %#ok<AGROW>
end
% Exact point-source manufacture on the real frozen dual graph; trace PV.
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
profile=a.state.dual_cfg.profile;src=[.2 -.1 -2];n=33;q=12;L=12;h=24/n;
[z,w]=deepsea_bie_gauss(q);c=-L+h/2+(0:n-1)*h;
x=reshape(c(:)+z.'*h/2,[],1);wx=repmat(w.'*h/2,n,1);wx=wx(:);
[xx,yy]=meshgrid(x,x);weight=reshape(wx*wx.',[],1);
geo=deepsea3d_tfe_geometry(xx,yy,profile);jac=sqrt(1+geo.ex(:).^2+geo.ey(:).^2);
ps=[xx(:),yy(:),geo.eta(:)];normal=[-geo.ex(:),-geo.ey(:),ones(numel(xx),1)]./jac;
pf=ps;pf(:,3)=0;nf=repmat([0 0 1],size(pf,1),1);
center=deepsea3d_tfe_geometry(0,0,profile);target=[0 0 center.eta];assert(target(3)>0);
[g,grad]=deepsea3d_point_field(ps,src,k);flux=sum(grad.*normal,2);
[gf,gradf]=deepsea3d_point_field(pf,src,k);
[S,D]=deepsea_bie_layers(target,ps,normal,weight.*jac,[flux,g],k,-1);
[Sf,Df]=deepsea_bie_layers(target,pf,nf,weight,[gradf(:,3),gf],k,-1);
unreplaced=D(2)-S(1)-(Df(2)-Sf(1));
inside=abs(ps(:,1))<h/2 & abs(ps(:,2))<h/2;
[S0,D0]=deepsea_bie_layers(target,ps(inside,:),normal(inside,:),weight(inside).*jac(inside),[flux(inside),g(inside)],k,-1);
[Sf0,Df0]=deepsea_bie_layers(target,pf(inside,:),nf(inside,:),weight(inside),[gradf(inside,3),gf(inside)],k,-1);
[duffy,wd]=deepsea_bie_duffy([0 0],h,[0 0],32);
gd=deepsea3d_tfe_geometry(duffy(:,1),duffy(:,2),profile);jd=sqrt(1+gd.ex.^2+gd.ey.^2);
pd=[duffy,gd.eta];nd=[-gd.ex,-gd.ey,ones(size(jd))]./jd;pd0=pd;pd0(:,3)=0;
[ud,du]=deepsea3d_point_field(pd,src,k);fd=sum(du.*nd,2);[ud0,du0]=deepsea3d_point_field(pd0,src,k);
[Sd,Dd]=deepsea_bie_layers(target,pd,nd,wd.*jd,[fd,ud],k,-1);
[Sdf,Ddf]=deepsea_bie_layers(target,pd0,repmat([0 0 1],size(pd0,1),1),wd,[du0(:,3),ud0],k,-1);
delta=unreplaced-(D0(2)-S0(1)-(Df0(2)-Sf0(1)))+(Dd(2)-Sd(1)-(Ddf(2)-Sdf(1)));
exact=deepsea3d_point_field(target,src,k);trace=exact+delta;
trace_error=abs(trace-.5*exact)/abs(exact);
% Geometry derivatives/area are independently checked through finite differences.
xy=[0 0;8.7 1.2;10.1 -2.3;11.7 .2];step=1e-5;
ge=deepsea3d_tfe_geometry(xy(:,1),xy(:,2),profile);
gx1=deepsea3d_tfe_geometry(xy(:,1)+step,xy(:,2),profile);gx0=deepsea3d_tfe_geometry(xy(:,1)-step,xy(:,2),profile);
gy1=deepsea3d_tfe_geometry(xy(:,1),xy(:,2)+step,profile);gy0=deepsea3d_tfe_geometry(xy(:,1),xy(:,2)-step,profile);
ex=(gx1.eta-gx0.eta)/(2*step);ey=(gy1.eta-gy0.eta)/(2*step);
area_fd=sqrt(1+ex.^2+ey.^2);area=sqrt(1+ge.ex.^2+ge.ey.^2);
geometry_error=max([abs(ex-ge.ex);abs(ey-ge.ey);abs(area_fd-area)]);
assert(geometry_error<1e-7);
output=struct('windows',windows,'curved_trace_relative',trace_error, ...
    'curved_trace_real',real(trace),'curved_trace_imag',imag(trace), ...
    'expected_trace_real',real(.5*exact),'expected_trace_imag',imag(.5*exact), ...
    'geometry_FD_error',geometry_error,'wall_s',toc(timer), ...
    'decision','STOP: physical reference gates not satisfied; no dual scattering or PM solve', ...
    'G1_target_PM','OPEN','WP2_formal',false);
save(file,'output','-v7');writetable(struct2table(windows),fullfile(outdir,'bie_window.csv'));
fid=fopen(fullfile(outdir,'bie_window_diagnosis.json'),'w');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
fprintf('BIE actual-graph manufactured PV error %.4g; geometry FD %.3g; wall%.3fs\n',trace_error,geometry_error,output.wall_s);
end
