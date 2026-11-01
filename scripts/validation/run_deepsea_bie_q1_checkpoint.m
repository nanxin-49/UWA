function receipt=run_deepsea_bie_q1_checkpoint(n,base_degree)
%RUN_DEEPSEA_BIE_Q1_CHECKPOINT One bounded 6 kHz planar manufactured solve.
% Validation-only Q1 nodal density; frozen 3-D half-space layer kernels.
arguments
    n (1,1) double {mustBeMember(n,[16 24])}
    base_degree (1,1) double {mustBeMember(base_degree,[8 12])}
end
assert(base_degree==8 || n==24,'bieQ1:Variant','Only n24 has a q12 variant.');
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outroot=fullfile(root,'results','validation','deepsea_bie_q1_20261026');
if ~isfolder(outroot),mkdir(outroot);end
out=fullfile(outroot,sprintf('n%d_q%d',n,base_degree));
assert(~isfolder(out),'bieQ1:Exists','No repeat or overwrite.');mkdir(out);
totalclock=tic;
k=2*pi*6000/1500;aux=-1;sigma=.7;L=4;
p=local_plan(L,n,base_degree,k,aux);
cp=struct('phase','assembled','n',n,'base_degree',base_degree,'L_m',L, ...
    'DOF',numel(p.interior_ids),'quadrature_sources',numel(p.weight), ...
    'configuration','Q1 continuous nodes; Gauss base, Duffy24 corner, Gauss16 neighbor; a=-1m', ...
    'elapsed_s',toc(totalclock));
[S,D]=deepsea_bie_layers(p.target,p.source,p.normal,p.weight,p.basis,k,aux);
Afull=D-1i*k*S;
clear S D;
nearclock=tic;
for row=1:numel(p.interior_ids)
    ix=p.target_ix(row);iy=p.target_iy(row);t=p.target(row,:);
    for cy=max(1,iy-2):min(n,iy+1)
        for cx=max(1,ix-2):min(n,ix+1)
            nodes=local_nodes(n,cx,cy);
            base=local_cell_integral(t,p,cx,cy,base_degree,0);
            if (cx==ix || cx==ix-1) && (cy==iy || cy==iy-1)
                refined=local_cell_integral(t,p,cx,cy,24,1);
            else
                refined=local_cell_integral(t,p,cx,cy,16,0);
            end
            Afull(row,nodes)=Afull(row,nodes)+refined-base;
        end
    end
end
near_wall_s=toc(nearclock);
A=Afull(:,p.interior_ids)+.5*eye(numel(p.interior_ids));
clear Afull;
cp.assembly_wall_s=toc(totalclock);cp.near_wall_s=near_wall_s;
save(fullfile(out,'assembled.mat'),'A','cp','-v7');

rhs=deepsea_bie_gaussian_reference(p.target,k,aux,sigma,512);
rhs256=deepsea_bie_gaussian_reference(p.target,k,aux,sigma,256);
reference_change=norm(rhs-rhs256)/norm(rhs);
solveclock=tic;phi=A\rhs;solve_wall_s=toc(solveclock);
fresh=norm(A*phi-rhs)/norm(rhs);
all_phi=complex(zeros((n+1)^2,1));all_phi(p.interior_ids)=phi;
node_truth=exp(-sum(p.node_xy.^2,2)/(2*sigma^2));
density_relative=norm(phi-node_truth(p.interior_ids))/norm(node_truth(p.interior_ids));
cp.phase='solved';cp.fresh_linear=fresh;cp.solve_wall_s=solve_wall_s;
cp.reference_change=reference_change;cp.elapsed_s=toc(totalclock);
save(fullfile(out,'solved.mat'),'phi','all_phi','rhs','cp','-v7');

shift=[.27 .19];
boundary_target=[p.target(:,1)+shift(1)*p.h,p.target(:,2)+shift(2)*p.h,zeros(numel(phi),1)];
boundary=local_boundary(p,all_phi,boundary_target,shift);
boundary_truth=deepsea_bie_gaussian_reference(boundary_target,k,aux,sigma,512);
boundary_relative=norm(boundary-boundary_truth)/norm(boundary_truth);
obs=[0 0 .3;.5 .2 1;0 0 5];
[field16,continuous16]=local_field(p,all_phi,obs,16,sigma);
[field24,~]=local_field(p,all_phi,obs,24,sigma);
field_truth=deepsea_bie_gaussian_reference(obs,k,aux,sigma,512);
field_relative=norm(field16-field_truth)/norm(field_truth);
field_integral_change=norm(field16-field24)/norm(field24);
known_density_field_error=norm(continuous16-field_truth)/norm(field_truth);
window5=local_continuous_field(5,round(5/p.h)*2,16,k,aux,sigma,obs);
known_window5_error=norm(window5-field_truth)/norm(field_truth);

central=ceil(n/2)+1;t=[p.x(central),p.x(central),0];
self24=local_cell_integral(t,p,central-1,central-1,24,1);
self32=local_cell_integral(t,p,central-1,central-1,32,1);
neighbor16=local_cell_integral(t,p,central+1,central-1,16,0);
neighbor24=local_cell_integral(t,p,central+1,central-1,24,0);
self_change=norm(self24-self32)/norm(self32);
neighbor_change=norm(neighbor16-neighbor24)/norm(neighbor24);

receipt=struct('frequency_hz',6000,'n',n,'h_m',p.h,'L_m',L,'DOF',numel(phi), ...
    'base_degree',base_degree,'self_Duffy_degree',24,'neighbor_Gauss_degree',16, ...
    'fresh_linear',fresh,'independent_boundary',boundary_relative, ...
    'independent_field',field_relative,'density_error',density_relative, ...
    'reference_change',reference_change,'self_Duffy24_32_change',self_change, ...
    'neighbor_Gauss16_24_change',neighbor_change, ...
    'field_Gauss16_24_change',field_integral_change, ...
    'known_gaussian_L4_field_error',known_density_field_error, ...
    'known_gaussian_L5_field_error',known_window5_error, ...
    'assembly_wall_s',cp.assembly_wall_s,'near_wall_s',near_wall_s, ...
    'solve_wall_s',solve_wall_s,'total_kernel_wall_s',toc(totalclock), ...
    'linear_pass',fresh<1e-8,'reference_pass',max(boundary_relative,field_relative)<1e-4, ...
    'G1_target_PM','OPEN','WP2_formal',false, ...
    'unknown_exterior_closure_verified',false);
save(fullfile(out,'result.mat'),'receipt','phi','all_phi','rhs','boundary','boundary_truth', ...
    'field16','field24','field_truth','node_truth','p','-v7');
fid=fopen(fullfile(out,'receipt.json'),'w');assert(fid>=0);holder=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt);
end

function p=local_plan(L,n,q,k,aux)
h=2*L/n;x=linspace(-L,L,n+1);
[xx,yy]=meshgrid(x,x);node_xy=[xx(:),yy(:)];
% Explicit ndgrid order keeps target and interior-ID order identical.
[IY,IX]=ndgrid(2:n,2:n);ix=IX(:);iy=IY(:);
ids=sub2ind([n+1,n+1],iy,ix);
target=[node_xy(ids,:),zeros(numel(ids),1)];
[v,w]=deepsea_bie_gauss(q);u=(v+1)/2;wu=w/2;
[U,V]=meshgrid(u,u);local=[U(:),V(:)];weight_local=reshape(wu*wu.',[],1)*h^2;
shape=local_shape(local);Q=size(local,1);M=n^2*Q;
source=zeros(M,3);weight=zeros(M,1);
rr=zeros(4*M,1);cc=rr;vv=rr;
for cx=1:n
    for cy=1:n
        cellno=sub2ind([n n],cy,cx);id=(cellno-1)*Q+(1:Q);
        source(id,1)=x(cx)+h*local(:,1);
        source(id,2)=x(cy)+h*local(:,2);
        weight(id)=weight_local;
        nodes=local_nodes(n,cx,cy);
        for b=1:4
            pos=(b-1)*M+id;rr(pos)=id;cc(pos)=nodes(b);vv(pos)=shape(:,b);
        end
    end
end
p=struct('L',L,'n',n,'h',h,'k',k,'aux',aux,'x',x,'node_xy',node_xy, ...
    'interior_ids',ids,'target',target,'target_ix',ix,'target_iy',iy, ...
    'source',source,'normal',repmat([0 0 1],M,1),'weight',weight, ...
    'basis',sparse(rr,cc,vv,M,(n+1)^2));
end

function nodes=local_nodes(n,cx,cy)
nodes=sub2ind([n+1,n+1],[cy cy cy+1 cy+1],[cx cx+1 cx cx+1]);
end

function shape=local_shape(uv)
u=uv(:,1);v=uv(:,2);
shape=[(1-u).*(1-v),u.*(1-v),(1-u).*v,u.*v];
end

function values=local_cell_integral(target,p,cx,cy,q,duffy)
h=p.h;left=[p.x(cx),p.x(cy)];center=left+h/2;
if duffy
    [xy,w]=deepsea_bie_duffy(center,h,target(1:2),q);
else
    [v,g]=deepsea_bie_gauss(q);u=(v+1)/2;
    [X,Y]=meshgrid(u,u);xy=[left(1)+h*X(:),left(2)+h*Y(:)];
    w=reshape(g*g.',[],1)*h^2/4;
end
shape=local_shape((xy-left)/h);
source=[xy,zeros(size(xy,1),1)];normal=repmat([0 0 1],size(xy,1),1);
[S,D]=deepsea_bie_layers(target,source,normal,w,shape,p.k,p.aux);
values=D-1i*p.k*S;
end

function boundary=local_boundary(p,all_phi,target,shift)
q=size(p.source,1)/p.n^2;
q=sqrt(q);
[S,D]=deepsea_bie_layers(target,p.source,p.normal,p.weight,p.basis*all_phi,p.k,p.aux);
boundary=D-1i*p.k*S;
for row=1:size(target,1)
    ix=p.target_ix(row);iy=p.target_iy(row);t=target(row,:);
    corner=local_nodes(p.n,ix,iy);
    interpolation=local_shape(shift);
    boundary(row)=boundary(row)+.5*(interpolation*all_phi(corner));
    for cy=max(1,iy-1):min(p.n,iy+1)
        for cx=max(1,ix-1):min(p.n,ix+1)
            nodes=local_nodes(p.n,cx,cy);
            base=local_cell_integral(t,p,cx,cy,q,0);
            if cx==ix && cy==iy
                refined=local_cell_integral(t,p,cx,cy,24,1);
            else
                refined=local_cell_integral(t,p,cx,cy,16,0);
            end
            boundary(row)=boundary(row)+(refined-base)*all_phi(nodes);
        end
    end
end
end

function [field,continuous]=local_field(p,phi,obs,q,sigma)
pp=local_plan(p.L,p.n,q,p.k,p.aux);
source_density=pp.basis*phi;
[S,D]=deepsea_bie_layers(obs,pp.source,pp.normal,pp.weight,source_density,p.k,p.aux);
field=D-1i*p.k*S;
continuous_density=exp(-sum(pp.source(:,1:2).^2,2)/(2*sigma^2));
[S,D]=deepsea_bie_layers(obs,pp.source,pp.normal,pp.weight,continuous_density,p.k,p.aux);
continuous=D-1i*p.k*S;
end

function field=local_continuous_field(L,n,q,k,aux,sigma,obs)
p=local_plan(L,n,q,k,aux);
density=exp(-sum(p.source(:,1:2).^2,2)/(2*sigma^2));
[S,D]=deepsea_bie_layers(obs,p.source,p.normal,p.weight,density,k,aux);
field=D-1i*k*S;
end
