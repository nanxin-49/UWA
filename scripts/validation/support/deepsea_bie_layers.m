function [S,D]=deepsea_bie_layers(target,source,normal,weight,density,k,aux,guard)
%DEEPSEA_BIE_LAYERS Streamed 3-D S and water-normal D, no dense global matrix.
% aux<min(surface) uses G-G(image); Inf denotes free-space diagnostic.
arguments
    target (:,3) double
    source (:,3) double
    normal (:,3) double
    weight (:,1) double
    density (:,:) double
    k (1,1) double
    aux (1,1) double
    guard = []
end
assert(size(source,1)==size(density,1) && size(source,1)==numel(weight));
S=complex(zeros(size(target,1),size(density,2)));D=S;
weighted=density.*weight;
for i=1:16:size(target,1)
    ti=i:min(i+15,size(target,1));
    for j=1:4096:size(source,1)
        if ~isempty(guard),guard();end
        sj=j:min(j+4095,size(source,1));
        [G,K]=kernel(target(ti,:),source(sj,:),normal(sj,:),k);
        if isfinite(aux)
            image=source(sj,:);image(:,3)=2*aux-image(:,3);
            reflected=normal(sj,:);reflected(:,3)=-reflected(:,3);
            [Gi,Ki]=kernel(target(ti,:),image,reflected,k);G=G-Gi;K=K-Ki;
        end
        S(ti,:)=S(ti,:)+G*weighted(sj,:);D(ti,:)=D(ti,:)+K*weighted(sj,:);
    end
end
end
function [G,K]=kernel(x,y,n,k)
dx=y(:,1).'-x(:,1);dy=y(:,2).'-x(:,2);dz=y(:,3).'-x(:,3);
R=sqrt(dx.^2+dy.^2+dz.^2);assert(all(R(:)>0),'BIE:Coincidence','Use self-panel Duffy quadrature.');
G=exp(1i*k*R)./(4*pi*R);
K=G.*(1i*k-1./R).*(dx.*n(:,1).'+dy.*n(:,2).'+dz.*n(:,3).')./R;
end
