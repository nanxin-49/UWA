function y=deepsea_bie_flat_action(phi,p,guard)
%DEEPSEA_BIE_FLAT_ACTION Corrected flat trace of D_a-i*k*S_a, water side.
if nargin<3,guard=[];end
[S,D]=deepsea_bie_layers(p.target,p.source,p.normal,p.weight,phi(p.cell_index),p.k,p.aux,guard);
y=.5*phi+D-1i*p.k*S+p.near_correction*phi;
end
