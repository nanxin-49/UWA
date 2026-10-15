function output=estimate_deepsea3d_bie_cost()
%ESTIMATE_DEEPSEA3D_BIE_COST Off-diagonal arithmetic only, no BIE solver.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
outdir=fullfile(root,'results','validation','deepsea_two_track_20261014');
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
cfg=a.state.pm_cfg;k=2*pi*cfg.frequency_hz/cfg.c0_mps;
x=linspace(-11,11,16);[xx,yy]=meshgrid(x,x);
g=deepsea3d_tfe_geometry(xx,yy,cfg.profile);
q=[xx(:),yy(:),g.eta(:)];normal=[-g.ex(:),-g.ey(:),ones(numel(xx),1)];
normal=normal./sqrt(sum(normal.^2,2));target=q+[.1 .07 .03];
times=zeros(1,6);checksum=0;
for repeat=1:6
    timer=tic;
    dx=target(:,1)-q(:,1).';dy=target(:,2)-q(:,2).';dz=target(:,3)-q(:,3).';
    R=sqrt(dx.^2+dy.^2+dz.^2);Green=exp(1i*k*R)./(4*pi*R);
    dn=-(dx.*normal(:,1).'+dy.*normal(:,2).'+dz.*normal(:,3).')./R;
    combined=Green.*((1i*k-1./R).*dn-1i*k);
    checksum=norm(combined*ones(size(q,1),1));times(repeat)=toc(timer);
end
pairs=size(q,1)^2;per_pair=median(times(2:end))/pairs;
rows=struct([]);
for h=[.25 .125 .0625 1/24]
    nodes=(round(24/h)+1)^2;dense=16*nodes^2/2^30;
    row=struct('patch_width_m',24,'h_m',h,'nodal_DOF',nodes, ...
        'dense_complex_matrix_GiB',dense,'pairs_per_direct_matvec',nodes^2, ...
        'block256_kernel_MiB',16*256*nodes/2^20, ...
        'arithmetic_only_matvec_s',nodes^2*per_pair);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
output=struct('scope','proposal only; not a BIE assembly, solve, or reference', ...
    'microkernel','off-diagonal free-space combined-layer arithmetic; no singular/near quadrature', ...
    'sample_pairs',pairs,'warm_times_s',times(2:end),'seconds_per_pair',per_pair, ...
    'checksum',checksum,'rows',rows,'G1_target_PM','OPEN','WP2_formal',false, ...
    'timing_limit','small-array arithmetic extrapolation; excludes quadrature, windows, exterior correction, iterations, overhead; not an affordability certificate');
writetable(struct2table(rows),fullfile(outdir,'BIE_cost_estimate.csv'));
fid=fopen(fullfile(outdir,'BIE_cost_estimate.json'),'w');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
fprintf('BIE arithmetic only %.3g ns/pair; no BIE solve performed\n',per_pair*1e9);
end
