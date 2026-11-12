function [u, grad] = deepsea3d_point_field(points, source, k)
%DEEPSEA3D_POINT_FIELD Unit Green monopole, exp(-i omega t), SI coordinates.
% (Delta+k^2)u=-delta; u=e^(ikR)/(4*pi*R). No calibrated source power.
arguments
    points (:,3) double {mustBeFinite}
    source (1,3) double {mustBeFinite}
    k (1,1) double {mustBePositive,mustBeFinite}
end
d=points-source; R=sqrt(sum(d.^2,2));
assert(all(R>0),'deepsea3d:SourceCoincidence','Evaluation at source.');
u=exp(1i*k*R)./(4*pi*R);
grad=(u.*(1i*k-1./R)./R).*d;
end
