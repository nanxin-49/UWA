function audit=audit_deepsea_bie_reconstruction()
%AUDIT_DEEPSEA_BIE_RECONSTRUCTION Saved-data, independent component/PDE checks.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
outdir=fullfile(root,'results','validation','deepsea_bie_reconstruction_20261014');
timer=tic;a=load(fullfile(outdir,'bie_minimum.mat'),'output','raw');
b=load(fullfile(outdir,'frequency_reconstruction.mat'),'output','saved');
d=load(fullfile(outdir,'bie_window_diagnosis.mat'),'output');
assert(all([a.output.cases.linear_pass]) && ~any([a.output.cases.reference_pass]));
assert(~a.output.PM_BIE_solve && ~a.output.reference_ready && ~a.output.WP2_formal);
assert(all([b.output.rows.passed]) && ~b.output.dual_intermediate_certified && ~b.output.generic_total_only_reconstruction);
components=zeros(numel(b.saved),1);
for j=1:numel(b.saved)
    s=b.saved(j);q=s.q;Ld=sqrt((1505*q)^2+1495^2);Lr=sqrt((1505*q)^2+1505^2);
    reflected=deepsea_delay_separated_reconstruct(s.source_f_hz,s.path_response(:,2),s.delay_s(2),s.query_hz,b.output.meta);
    exact=-exp(-2i*pi*s.query_hz*s.delay_s(2))/(4*pi*Lr);
    components(j)=norm(reflected-exact)/norm(exact);assert(components(j)<1e-3);
    assert(norm(s.reconstructed-s.exact)/(sqrt(numel(s.query_hz))/(4*pi*Ld))<1e-3);
end
% Independent finite differences of the actually solved layer field, not its RHS.
phi=a.raw{end}.phi;k=2*pi*6000/1500;p=deepsea_bie_flat_plan(3,40,k,-1);
origin=[.3 .2 1];base=field(origin);pde=zeros(1,2);
for j=1:2
    step=(1500/6000)/256/2^(j-1);lap=0;
    for axis=1:3
        off=zeros(1,3);off(axis)=step;
        lap=lap+(field(origin+off)-2*base+field(origin-off))/step^2;
    end
    pde(j)=abs(lap+k^2*base)/(k^2*abs(base));
end
ray=struct([]);step=(1500/6000)/512;
for r=[20 40 80]
    center=field([0 0 r]);derivative=(field([0 0 r+step])-field([0 0 r-step]))/(2*step);
    row=struct('r_m',r,'outgoing_mismatch',abs(derivative-1i*k*center)/(k*abs(center)), ...
        'incoming_mismatch',abs(derivative+1i*k*center)/(k*abs(center)));
    if isempty(ray),ray=row;else,ray(end+1)=row;end %#ok<AGROW>
end
assert(max(pde)<1e-4 && all([ray.incoming_mismatch]>1.5));
audit=struct('status','AUDIT_PASS_WITH_REFERENCE_REJECTION','all_six_linear_pass',true, ...
    'all_six_reference_fail',true,'reflected_component_errors',components, ...
    'layer_field_PDE_FD',pde,'outgoing_ray',ray, ...
    'radiation_scope','one-ray FD consistency plus outgoing kernels; not a complete exterior-domain certificate', ...
    'A_total_wall_including_local_checks_s',a.output.total_wall_s+d.output.wall_s, ...
    'B_wall_s',b.output.wall_s,'audit_wall_s',toc(timer), ...
    'A_decision','STOP_REFERENCE_GATES_FAIL; no dual scattering / PM solve', ...
    'B_decision','PASS_ANALYTIC_KNOWN_PATHS_ONLY','generic_frequency_reconstruction_accepted',false, ...
    'G1_target_PM','OPEN','WP2_formal',false);
fid=fopen(fullfile(outdir,'data_audit.json'),'w');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
fprintf('Data audit: reflected %.3g; PDE FD %s; A rejected, analytic B pass\n',max(components),mat2str(pde,4));
    function value=field(points)
        [S,D]=deepsea_bie_layers(points,p.source,p.normal,p.weight,phi(p.cell_index),k,-1);
        value=D-1i*k*S;
    end
end
