function p=deepsea_bie_flat_plan(L,n,k,aux,offset)
%DEEPSEA_BIE_FLAT_PLAN Piecewise-constant cells, corrected singular/near cells.
if nargin<5,offset=[0 0];end
h=2*L/n;x=-L+h/2+(0:n-1)*h;[xx,yy]=meshgrid(x,x);center=[xx(:),yy(:)];
[z,w]=deepsea_bie_gauss(4);[dx,dy]=meshgrid(z*h/2,z*h/2);ww=(w*w.')*h^2/4;
count=n^2;q=numel(ww);
sx=center(:,1)+dx(:).';sy=center(:,2)+dy(:).';
source=[sx(:),sy(:),zeros(count*q,1)];normal=repmat([0 0 1],size(source,1),1);
weight=repelem(ww(:),count);cell_index=repmat((1:count).',q,1);
target=[center+h*offset,zeros(count,1)];ri=[];ci=[];val=[];
for ix=-1:1
    for iy=-1:1
        c=h*[ix iy];txy=h*offset;
        if ix==0 && iy==0
            [near,wn]=deepsea_bie_duffy(c,h,txy,24);
        else
            [v,wv]=deepsea_bie_gauss(16);[a,b]=meshgrid(v*h/2,v*h/2);
            near=[c(1)+a(:),c(2)+b(:)];wn=reshape(wv*wv.'*h^2/4,[],1);
        end
        sn=[near,zeros(size(near,1),1)];nn=repmat([0 0 1],size(sn,1),1);
        [S,D]=deepsea_bie_layers([txy 0],sn,nn,wn,ones(size(wn)),k,aux);
        sb=[c(1)+dx(:),c(2)+dy(:),zeros(q,1)];
        [S0,D0]=deepsea_bie_layers([txy 0],sb,repmat([0 0 1],q,1),ww(:),ones(q,1),k,aux);
        correction=D-1i*k*S-(D0-1i*k*S0);
        [row,col]=ndgrid(1:n,1:n);inside=row+iy>=1 & row+iy<=n & col+ix>=1 & col+ix<=n;
        from=sub2ind([n n],row(inside),col(inside));
        to=sub2ind([n n],row(inside)+iy,col(inside)+ix);
        ri=[ri;from];ci=[ci;to];val=[val;repmat(correction,numel(from),1)]; %#ok<AGROW>
    end
end
p=struct('L',L,'n',n,'h',h,'k',k,'aux',aux,'center',center,'target',target, ...
    'source',source,'normal',normal,'weight',weight,'cell_index',cell_index, ...
    'near_correction',sparse(ri,ci,val,count,count));
end
