function tests=test_deepsea_bie_reconstruction
tests=functiontests(localfunctions);
end
function setupOnce(~)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_vertical_project;
end
function testGaussIndependentMoments(t)
[x,w]=deepsea_bie_gauss(6);
for power=0:11
    expected=0;if mod(power,2)==0,expected=2/(power+1);end
    verifyLessThan(t,abs(sum(w.*x.^power)-expected),1e-13);
end
end
function testLayerReciprocityAndBoundaryImageZero(t)
x=[.2 .7 1.3];y=[-.4 .1 .2];normal=[0 0 1];k=3.2;
[a,~]=deepsea_bie_layers(x,y,normal,1,1,k,-1);
[b,~]=deepsea_bie_layers(y,x,normal,1,1,k,-1);
verifyLessThan(t,abs(a-b),1e-15);
[zero,~]=deepsea_bie_layers([.2 .7 -1],y,normal,1,1,k,-1);
verifyLessThan(t,abs(zero),1e-15);
end
function testLaplaceSelfHasExactIndependentSquareIntegral(t)
[xy,w]=deepsea_bie_duffy([0 0],1,[0 0],24);
[S,~]=deepsea_bie_layers([0 0 0],[xy zeros(size(xy,1),1)], ...
    repmat([0 0 1],size(xy,1),1),w,ones(size(w)),0,Inf);
verifyLessThan(t,abs(S-asinh(1)/pi),1e-12);
end
function testCoincidentKernelFailsInsteadOfDroppingSelf(t)
verifyError(t,@()deepsea_bie_layers([0 0 0],[0 0 0],[0 0 1],1,1,2,-1),'BIE:Coincidence');
end
function testDelayReconstructionCannotCertifyDual(t)
meta=struct('physics_tier',"compact_dual_conditional",'known_path_components',true, ...
    'G1_target_PM',"OPEN",'WP2_formal',false);
verifyError(t,@()deepsea_delay_separated_reconstruct([4000;8000],ones(2,1),0,6000,meta), ...
    'Reconstruction:Scope');
end
function testKnownDelayRecoversIndependentOffGridExponential(t)
f=(4000:1000:8000).';q=(4123:713:7999).';tau=.006666666666667;
meta=struct('physics_tier',"analytic_control",'known_path_components',true, ...
    'G1_target_PM',"OPEN",'WP2_formal',false);
sample=exp(-2i*pi*f*tau);actual=deepsea_delay_separated_reconstruct(f,sample,tau,q,meta);
verifyLessThan(t,norm(actual-exp(-2i*pi*q*tau)),1e-12);
end
