function receipt = audit_bellhop_rotated_field_primitives(transport_dir,chart_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_FIELD_PRIMITIVES Replay exact Q/SHD/coordinate evidence.
arguments
    transport_dir (1,:) char
    chart_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
diary(fullfile(output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
bindings=struct('path',{},'sha256',{});
a=load(fullfile(transport_dir,'result.mat'),'runs');
b=load(fullfile(chart_dir,'result.mat'),'runs');
all_runs=[a.runs b.runs];
for j=1:numel(all_runs)
    run=all_runs{j};
    if j<=5
        names={'uniform_native','uniform_point','munk_linear_native','munk_fast_native','munk_point'};
        folder=fullfile(transport_dir,names{j});
    else
        names={'uniform_right','munk_right_native','munk_right_point','munk_left_point','munk_up_point','munk_down_point'};
        folder=fullfile(chart_dir,names{j-5});
    end
    stem=fullfile(folder,'case'); raw=read_bellhop_rotated_field(stem);
    add(sprintf('case%d_raw_SHD_MAT_pressure',j),max(abs(raw.pressure-run.field.pressure),[],'all'),0);
    add(sprintf('case%d_exact_range_binding',j),max(abs(raw.receiver_range_m-run.cfg.receiver_range_m(:))),0);
    for ext={'.env','.ssp','.iwcov','.r25d','.shd','.r25diag','.r25sample','.prt'}
        path=[stem ext{1}]; bindings(end+1)=struct('path',path,'sha256',black_file_hash(path));
    end
end
% Direct-axis ray-theory reference for the exact written linear Q medium.
run=a.runs{5}; folder=fullfile(transport_dir,'munk_point');
fid=fopen(fullfile(folder,'case.ssp'),'r'); nr=fscanf(fid,'%d',1);
rr=1000*fscanf(fid,'%f',nr); cc=fscanf(fid,'%f',nr); fclose(fid);
R=run.cfg.receiver_range_m; amplitude=zeros(size(R)); tau=zeros(size(R));
for j=1:numel(R)
    nodes=[rr(rr<R(j));R(j)]; c=interp1(rr,cc,nodes);
    I=trapz(nodes,c); amplitude(j)=sqrt(cc(1)*c(end))/I;
    u=diff(c)./c(1:end-1); ratio=ones(size(u)); nonzero=u~=0;
    ratio(nonzero)=log1p(u(nonzero))./u(nonzero);
    tau(j)=sum(diff(nodes)./c(1:end-1).*ratio);
end
pressure=-conj(run.field.pressure);
add('exact_Q_3D_transport_amplitude',max(abs(abs(pressure)./amplitude-1)),2e-6);
add('exact_Q_ray_phase_rad',max(abs(angle(pressure.*exp(-1i*2*pi*10000*tau)))),1e-4);
uniform=-conj(a.runs{2}.field.pressure); r=a.runs{2}.cfg.receiver_range_m;
expected=exp(1i*2*pi*10000*r/1500)./r;
add('exact_coordinate_uniform_Green',norm(uniform-expected)/norm(expected),1e-6);
metadata_offset=max(abs(a.runs{2}.field.receiver_range_m-a.runs{2}.field.receiver_range_metadata_m));
metadata_phase_bound=2*pi*10000*metadata_offset/1500;
% Source plane is independently computed by both horizontal charts.
Q=b.runs{1}.cfg.receiver_q_m; R=b.runs{1}.cfg.receiver_range_m;
distance=hypot(Q(:),R); exact=exp(1i*2*pi*10000*distance/1500)./distance;
p=-conj(b.runs{1}.field.pressure); valid=b.runs{1}.coverage;
add('source_chart_uniform_Green',norm(p(valid)-exact(valid))/norm(exact(valid)),1e-3);
center=find(Q==0);
add('source_plane_non_source_points_finite',double(any(~isfinite([b.runs{3}.field.pressure(center,:),b.runs{4}.field.pressure(center,:)]))),0);
for j=1:4
    r=R(j); ip=find(abs(Q-r)<1e-12); im=find(abs(Q+r)<1e-12);
    values=[b.runs{5}.field.pressure(ip,j),b.runs{3}.field.pressure(im,j); ...
        b.runs{5}.field.pressure(im,j),b.runs{4}.field.pressure(ip,j); ...
        b.runs{6}.field.pressure(im,j),b.runs{3}.field.pressure(ip,j); ...
        b.runs{6}.field.pressure(ip,j),b.runs{4}.field.pressure(im,j)];
    add(sprintf('chart_overlap_r%.2g_relative',r),max(abs(values(:,1)-values(:,2))./abs(values(:,2))),1e-3);
end
receipt=struct('status','COMPLETE','failed_checks',nnz(~checks.pass),'elapsed_s',toc(t), ...
    'scope','Raw coordinate/pressure binding, exact written-Q ray theory and direct chart overlap; not independent FFP/full field', ...
    'native_metadata_max_offset_m',metadata_offset,'native_metadata_phase_error_bound_rad',metadata_phase_bound, ...
    'matlab_release',version,'result_dir',output_dir);
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'checks','receipt','bindings','amplitude','tau');
fid=fopen(fullfile(output_dir,'input_sha256.json'),'w'); fprintf(fid,'%s',jsonencode(bindings,PrettyPrint=true)); fclose(fid);
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
