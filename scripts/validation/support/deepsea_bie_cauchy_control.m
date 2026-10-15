function rows=deepsea_bie_cauchy_control(profile)
%DEEPSEA_BIE_CAUCHY_CONTROL Exact infinite-flat background + compact difference.
% Known point-source Cauchy data test, not an unknown-density scattering solve.
k=2*pi*6000/1500;aux=-1;src=[.2 -.1 -2];rows=struct([]);
targets=[0 0 1;1 .2 2;0 0 5];exact=deepsea3d_point_field(targets,src,k);
for n=[17 25 33]
    for q=[8 12]
        timer=tic;L=12;h=2*L/n;centers=-L+h/2+(0:n-1)*h;
        [z,w]=deepsea_bie_gauss(q);x=reshape(centers(:)+(z.'*h/2),[],1);
        wx=repmat(w.'*h/2,n,1);wx=wx(:);[xx,yy]=meshgrid(x,x);weight=reshape(wx*wx.',[],1);
        geo=deepsea3d_tfe_geometry(xx,yy,profile);
        jac=sqrt(1+geo.ex(:).^2+geo.ey(:).^2);
        normal=[-geo.ex(:),-geo.ey(:),ones(numel(xx),1)]./jac;
        points=[xx(:),yy(:),geo.eta(:)];flat=points;flat(:,3)=0;
        [g,grad]=deepsea3d_point_field(points,src,k);flux=sum(grad.*normal,2);
        [gf,gradf]=deepsea3d_point_field(flat,src,k);nf=repmat([0 0 1],size(flat,1),1);
        [S,D]=deepsea_bie_layers(targets,points,normal,weight.*jac,[flux,g],k,aux);
        [Sf,Df]=deepsea_bie_layers(targets,flat,nf,weight,[gradf(:,3),gf],k,aux);
        delta=D(:,2)-S(:,1)-(Df(:,2)-Sf(:,1));
        full=exact+delta;raw=D(:,2)-S(:,1);
        row=struct('panels_per_side',n,'gauss_degree',q,'source_points',numel(weight), ...
            'compact_difference_relative',norm(delta)/norm(exact), ...
            'raw_patch_relative',norm(raw-exact)/norm(exact), ...
            'corrected_relative',norm(full-exact)/norm(exact), ...
            'min_eta',min(geo.eta(:)),'max_eta',max(geo.eta(:)), ...
            'normal_unit_error',max(abs(sum(normal.^2,2)-1)),'wall_s',toc(timer));
        if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    end
end
end
