function tests=test_deepsea3d_controls
tests=functiontests(localfunctions);
end
function testPointNormalizationAndDerivative(t)
k=2*pi*6000/1500; src=[0 0 1500]; p=[100 200 0;300 400 5];
[u,d]=deepsea3d_point_field(p,src,k); R=vecnorm(p-src,2,2);
verifyLessThan(t,max(abs(abs(u).*4*pi.*R-1)),1e-14);
h=1e-5;
for j=1:3
    a=p; b=p; a(:,j)=a(:,j)+h; b(:,j)=b(:,j)-h;
    fd=(deepsea3d_point_field(a,src,k)-deepsea3d_point_field(b,src,k))/(2*h);
    verifyLessThan(t,norm(fd-d(:,j))/norm(d(:,j)),2e-6);
end
end
function testSurfaceNormals(t)
p=local_profile; x=[0 5 9 11 13]; y=[0 2 3 1 4]; h=1e-5;
[~,ex,ey]=deepsea3d_surface(x,y,p);
fx=(deepsea3d_surface(x+h,y,p)-deepsea3d_surface(x-h,y,p))/(2*h);
fy=(deepsea3d_surface(x,y+h,p)-deepsea3d_surface(x,y-h,p))/(2*h);
verifyLessThan(t,max(abs(fx-ex)),1e-10); verifyLessThan(t,max(abs(fy-ey)),1e-10);
end
function testFlatBackground(t)
c=local_cfg; c.profile.amplitudes_m=0;
u=solve_deepsea3d_fullfk(c); b=deepsea3d_flat(c,0);
verifyEqual(t,u.reflect_native,b.reflect_native,'AbsTol',1e-15);
verifyEqual(t,u.rough_correction,0,'AbsTol',1e-15);
end
function testShiftedPlane(t)
c=local_cfg; c.profile.amplitudes_m=0; c.profile.height_m=.2;
u=solve_deepsea3d_fullfk(c); b=deepsea3d_flat(c,.2);
verifyEqual(t,u.reflect_native,b.reflect_native,'AbsTol',1e-15);
verifyEqual(t,b.excess_delay_s,(10-.4)/1500,'AbsTol',1e-15);
end
function testRotation(t)
c=local_cfg; c.profile=local_profile;
a=solve_deepsea3d_fullfk(c);
c.profile.wavevectors_radpm=[-c.profile.wavevectors_radpm(:,2) c.profile.wavevectors_radpm(:,1)];
b=solve_deepsea3d_fullfk(c);
verifyLessThan(t,abs(a.reflect_native-b.reflect_native)/abs(a.reflect_native),1e-11);
end
function testShapeReferenceZero(t)
c=local_cfg; c.profile.amplitudes_m=0;
r=solve_deepsea3d_shape_reference(c); verifyEqual(t,r.rough_correction,0,'AbsTol',1e-15);
end
function testQuadratureIdentity(t)
c=local_cfg; u=solve_deepsea3d_fullfk(c);
verifyLessThan(t,u.quadrature_identity_residual/abs(u.flat_exact),1e-12);
verifyLessThan(t,abs(u.rough_correction-u.term1_correction-u.term2_correction)/abs(u.flat_exact),1e-12);
end
function testPhaseBridge(t)
c=local_cfg; b=deepsea3d_flat(c,0);
verifyEqual(t,b.total_dsp,b.direct_dsp+b.reflect_dsp,'AbsTol',1e-15);
verifyEqual(t,b.direct_dsp,1/(4*pi*b.Ld_m),'AbsTol',1e-14);
verifyEqual(t,b.reflect_dsp,-exp(-1i*2*pi*c.frequency_hz*b.excess_delay_s)/(4*pi*b.Lr_m),'AbsTol',1e-14);
end
function c=local_cfg
c=struct('frequency_hz',6000,'c0_mps',1500,'tx_m',[0 0 1500], ...
    'rx_m',[0 0 5],'center_m',[0 0],'half_width_m',16,'spacing_m',.1, ...
    'quadrature','trap','profile',local_profile);
end
function p=local_profile
p=struct('height_m',0,'plateau_m',8,'support_m',12,'amplitudes_m',[.002 .001], ...
    'wavevectors_radpm',[.4 0;.2 .6],'phases_rad',[.2 .7]);
end
