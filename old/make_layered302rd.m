% make profiles across the domain.
% range-dependent layered


y1=-6000:100:6000; %transverse dimension
z1= -0:5:400;
x=0:20:30000;  % range-stepping dimension

[z,y]=meshgrid(z1,y1);

%cw= 1500 - sqrt((Pis*350/300));
wid=8; % layer width factor,smaller is sharper interface
zla=20; % layer center
cw1= 1500 + 12* erfc((z/wid)-(zla/wid))-(12*2) ;

cw=zeros( numel(x), size(cw1,1), size(cw1,2) ); % dim nx ny nz
n=1:numel(x);

zout=zeros( numel(x), size(z,1), size(z,2) );

m=ones(size(z(:,3:8)));
for ii=1:size(m,2)
    m(:,ii)=m(:,ii)*sech(ii-3).^1;
end

% to speed up gridded interpolant in the Envfile, need to keep the y,z grid
% fixed here. So alter cw instead.

for ii=n 
     cw(ii,:,:)=cw1;
     zout(ii,:,:)=z;
     %zout(ii,:,3:8)= z(:,3:8) +m.*2*sin(.8*ii/pi) ; 
     cw(ii,:,3:8)= cw1(:,3:8)-m.*3*sin(.8*ii/pi).^2; % some strange adhoc wiggles in x,z;
end

contourf( repmat(x,81,1).', squeeze(zout(:,10,:)),squeeze(cw(:,10,:)),[ min(cw(:)):.5:max(cw(:))] );
colorbar
%pcolor( y,z,squeeze(cw(1,:,:)))

save layered302rd  x y y1 zout cw


