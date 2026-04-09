% adapted from envF.m

% Y and Z make up the computation grid;
% Z is symmetric to enable image source
%
% Makes cin, mostly, from which U is made, then phase screen.  U=(cin-c0)./cin; 

eta_pert_ampl=2;  % Amplitude of perturbation of depths where cw is specified
                  % If zero, sound speed cw is as specified in filen.mat

if(~exist('cwi'))  % read in the slice, initially

dat=load(filen);  % z,y slices of sound speed, 
end


%iny=dat.x-mean(dat.x(:,1)); for test301 and test267, those have no x < 0.

%randomize the layering depths, quick internal wave perturbations.
%  pert=0*dat.z;
%  pert= eta_pert_ampl * randn(size(pert));
%  pert(:,1:2)=0 ; % sea surface
%  dat.z=dat.z+pert;

gridystep=Y(15,4)-Y(15,3);
ystretch=1;

% if(abs(waveorigin*ystretch) > max(max(Y)) )
%     yin=[ (waveorigin*ystretch):gridystep:(-waveorigin*ystretch)];  %vw
% else
%     yin=[(min(min(y))*ystretch):gridystep:(max(max(y))*ystretch) ];
% end

yin=[(min(min(y))*ystretch):gridystep:(max(max(y))*ystretch) ];

nyw=length(yin); Y2=ones(nz,1)*yin;
zgr=Z(:,1);
zin=[ min(Z(:,1))*1.1 ; zgr ; max(Z(:,1))*1.1]';

bthick =(max(Z(:))-depP)*1.1;
li=ones(size(yin));
bo=zin(min(find(zin > depP)));
wlast=zin(max(find(zin < depP)));

cb1=1650;  bcgr=1;   
cb2= cb1+bcgr*bthick;

      sea=find( Z(:,1) >=0 );
         
      % for range dependent, find the correct slices, interpolate
       
     if(dista ==0)
          
         cwslice=squeeze(dat.cw(1,:,:));
         zslice=squeeze(dat.zout(1,:,:)); 
         tt=cputime;
         cwi=griddata(dat.y,zslice,cwslice,Y(sea,:),Z(sea,:)); % (x z c xi zi)
         F=griddedInterpolant(dat.y,zslice,cwslice);
         tt=cputime; 
         cwi2=F(Y(sea,:).',Z(sea,:).').'; % (x z c xi zi)
         
         disp( ['cputime to grid cw at ' num2str(dista) ' = ' num2str(cputime-tt)])
     else
          
       a=find( dat.x-dista >=0 ,1) ; % find the slice past where we are
       
       if(dista > 0 &rdcw_flag ==1)
       if(a >1)
           d1=abs(dista-dat.x(a-1)) ; d2=abs(dista-dat.x(a)); 
       cwslice=(squeeze(dat.cw(a-1,:,:) ) *d1 + squeeze(dat.cw(a-1,:,:) ) *d2 )/(d1+d2);
       zslice=(squeeze(dat.zout(a-1,:,:) ) *d1 + squeeze(dat.zout(a-1,:,:) ) *d2 )/(d1+d2);
       else
           cwslice=squeeze(dat.cw(1,:,:));
           zslice=squeeze(dat.zout(1,:,:)); 
       end
        %tt=cputime;
    %          figure(100)
    %          imagesc(cwslice-F.Values)
    %          title(['cwnew-cwold, slice ' num2str(a)])
    %          colorbar
    %          pause(2)
    %          figure(1)
         F.Values=( cwslice ) ;
         cwi=F(Y(sea,:).',Z(sea,:).').';
         %cwi=griddata(cwslice,Y(sea,:),Z(sea,:)); % (x z c xi zi)
         %disp( ['cputime to grid cw at ' num2str(dista) ' = ' num2str(cputime-tt)])
         
       end
       
        end
  
    
    % domain is 15 percent flat
    
    dep= [ zeros(1, round(ny*.15))  ( Y(1,:)-min(Y(1,:) ))*100/max(Y(:))   ];
    dep=dep(1:ny);
   
    dz= abs(Z(3,5)-Z(4,5));
    dc=bcgr*dz;
    
    for ii=1:ny
        iz1=find( Z(sea,1) > dep(ii) ,1); 
        il=numel(cwi(iz1:end,ii) );
        cwi(iz1:end,ii)=cb1+linspace(0,il*dc,il);  
    end 
    
   
    
    cin= [flipud(cwi) ; cwi] ;            %  add in the image 
      
   
%bot= find(abs(Z) > dep);
bot=find(cin > 1600);
loss= (-.1) ; % db/lambda
fracperlamb= 10^(loss/20);
attb=fracperlamb.^(steplength/lambda);
attbo=ones(size(Z));
attbo(bot)=attb;

return

