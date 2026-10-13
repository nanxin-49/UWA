function W=deepsea3d_mapped_l0_solve(F,bottom,top,p)
%DEEPSEA3D_MAPPED_L0_SOLVE Exact discrete planar Robin inverse, no Taylor.
n=p.n;M=numel(p.s);
F_f=reshape(fft2(reshape(F,n,n,M)),n^2,M).';
g=fft2(deepsea3d_fg_resize(bottom,n));h=fft2(deepsea3d_fg_resize(top,n));
V_f=deepsea3d_tfe_modal_solve(p.modal,F_f,g(:),h(:));
W=reshape(ifft2(reshape(V_f.',n,n,M)),n^2,M);
end
