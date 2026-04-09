% make profiles across the domain.
% range-dependent layered


y1=-6000:100:6000; %transverse dimension
z1= -0:5:400;
x=[ 0 30000];  % range-stepping dimension

[z,y]=meshgrid(z1,y1);

%cw= 1500 - sqrt((Pis*350/300));
wid=8; % layer width factor,smaller is sharper interface
zla=20; % layer center
cw1= 1500 + 12* erfc((z/wid)-(zla/wid))-(12*2) ;

cw=zeros( numel(x), size(cw1,1), size(cw1,2) ); % dim nx ny nz
n=1:numel(x);

zout=zeros( numel(x), size(z,1), size(z,2) );

m=ones(size(z(:,2:8)));

for ii=n 
     cw(ii,:,:)=cw1;
     zout(ii,:,:)=z;
    % zout(ii,:,2:8)= z(:,2:8) +m.*ampl*sin(.8*ii/pi) ; 
end

%pcolor( y,z,squeeze(cw(1,:,:)))

save layered302ri  x y y1 zout cw


