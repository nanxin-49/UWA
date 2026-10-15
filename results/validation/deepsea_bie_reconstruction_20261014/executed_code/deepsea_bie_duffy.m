function [xy,weight]=deepsea_bie_duffy(center,h,target_xy,degree)
%DEEPSEA_BIE_DUFFY Four triangles from the target projection to a square.
[u,w]=deepsea_bie_gauss(degree);u=(u+1)/2;w=w/2;[U,V]=ndgrid(u,u);W=w*w.';
corners=center+h/2*[-1 -1;1 -1;1 1;-1 1];xy=[];weight=[];
for j=1:4
    a=corners(j,:)-target_xy;b=corners(mod(j,4)+1,:)-target_xy;
    ray=(1-V(:)).*a+V(:).*b;
    xy=[xy;target_xy+U(:).*ray]; %#ok<AGROW>
    weight=[weight;W(:).*U(:)*abs(det([a;b]))]; %#ok<AGROW>
end
end
