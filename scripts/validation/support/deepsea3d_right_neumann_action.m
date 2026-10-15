function y=deepsea3d_right_neumann_action(v,action)
%DEEPSEA3D_RIGHT_NEUMANN_ACTION Fixed right preconditioner P=2I-B.
% ACTION is the frozen L0-preconditioned B, never a modified physical PDE.
y=2*v-action(v);
end
