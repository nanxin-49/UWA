function top=deepsea3d_fg_top_product(Wtop,p,kind)
%DEEPSEA3D_FG_TOP_PRODUCT Includes the eta*T(Vprev) DtN product in dealiasing.
T=ifft2(1i*p.gamma.*fft2(Wtop));
if string(kind)=="galerkin"
    top=deepsea3d_fg_resize(-p.pad_geo.eta/p.b.*deepsea3d_fg_resize(T,2*p.n),p.n);
else
    top=deepsea3d_fg_resize(-p.geo.eta/p.b.*T,p.n);
end
end
