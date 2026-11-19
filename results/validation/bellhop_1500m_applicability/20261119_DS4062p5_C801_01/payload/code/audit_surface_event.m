function geom=audit_surface_event(root,run_dir)
% Independent flat-wall mirror geometry and logged physical pressure change.

iw=readlog([root '.iwdiag']);mw=readlog([root '.mwdiag']);d=readlog([root '.r25diag']);assert(height(iw)==801&&height(mw)==801&&height(d)==801);
geom=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
a=deg2rad(iw.alpha_deg);hitExpected=1490*tan(a);inc=[cos(a),sin(a)];refl=[-cos(a),sin(a)];rot=[cos(a),-sin(a)];
gate('all_first_sea_counts',nnz(mw.requested_order~=1|mw.walls_done~=1|mw.surface_bounces~=1|mw.bottom_bounces~=0),0);
gate('terminal_physical_bounces',nnz(d.top_bounces~=1|d.bottom_bounces~=0),0);
gate('reflection_wall_physical_depth_m',max(abs(iw.hit_r-1490)),1e-4);
gate('reflection_point_x_m',max(abs(iw.hit_z-hitExpected)),1e-4);
gate('incident_direction_rad',max(direrr([iw.inc_ur iw.inc_uz],inc)),1e-6);gate('reflected_direction_rad',max(direrr([iw.ref_ur iw.ref_uz],refl)),1e-6);gate('proper_rotated_direction_rad',max(direrr([iw.rot_ur iw.rot_uz],rot)),1e-6);
gate('surface_path_length_m',max(abs(1500*iw.tau_receiver_real-1496./cos(a))),1e-4);gate('surface_arrival_range_time_s',max(abs(iw.tau_receiver_real-1496./cos(a)/1500)),1e-7);gate('time_to_wall_s',max(abs(iw.tau_wall_real-1490./cos(a)/1500)),1e-7);
coeff=(iw.amp_ref./iw.amp_in).*exp(1i*(iw.phase_ref-iw.phase_in));gate('pressure_release_R_minus1',max(abs(coeff+1)),1e-12);gate('vacuum_amplitude_preserved',max(abs(iw.amp_ref-iw.amp_in)),1e-12);gate('pressure_release_phase_pi_rad',max(abs(angle(exp(1i*(iw.phase_ref-iw.phase_in-pi))))),1e-12);
gate('flat_surface_curvature',max(abs(iw.kappa)),0);gate('segment_membership',nnz(iw.wall_seg~=1&iw.wall_seg~=2),0);gate('segment_parameter_range',max([0;-iw.wall_lambda;iw.wall_lambda-1]),1e-10);
seg=round(iw.wall_seg);nodes=[-1000;0;1000];qlambda=nodes(seg)+1000*iw.wall_lambda;gate('segment_parameter_point_m',max(abs(qlambda-iw.hit_z)),1e-4);
writetable(iw,fullfile(run_dir,'all801_first_surface_states.csv'));writetable(geom,fullfile(run_dir,'surface_geometry_checks.csv'));
 function gate(name,value,limit),geom(end+1,:)={string(name),isscalar(value)&&isfinite(value)&&value<=limit,value,limit};end
end
function t=readlog(file)
fid=fopen(file,'r');
assert(fid>=0);
closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
hdr=fgetl(fid);
names=strsplit(strtrim(regexprep(hdr,'^#','')));
data=fscanf(fid,'%f',[numel(names) inf]).';
t=array2table(data,'VariableNames',names);
end
function e=direrr(u,v)
e=atan2(abs(u(:,1).*v(:,2)-u(:,2).*v(:,1)),sum(u.*v,2));
end
