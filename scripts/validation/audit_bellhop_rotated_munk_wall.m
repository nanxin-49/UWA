function receipt = audit_bellhop_rotated_munk_wall(input_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_MUNK_WALL Raw outputs vs exact written-Q Snell integrals.
arguments
    input_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
diary(fullfile(output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
data=load(fullfile(input_dir,'result.mat'),'results','receipt'); assert(strcmp(data.receipt.status,'COMPLETE'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
bindings=struct('path',{},'sha256',{}); reconstructed=cell(1,2);
for id=1:2
    run=data.results{id}; cfg=run.cfg; tag=cfg.reflection_wall;
    p=cell(1,2);
    for ib=1:2
        folder=fullfile(input_dir,sprintf('%s_mode%d',tag,ib));
        for ext={'.env','.ssp','.iwcov','.rotq','.shd','.iwdiag','.prt'}
            path=fullfile(folder,['case' ext{1}]);
            bindings(end+1)=struct('path',path,'sha256',black_file_hash(path));
        end
        p{ib}=read_bellhop_shd_unfolded_vertical(fullfile(folder,'case.shd'));
        add(sprintf('%s_raw_SHD_MAT_mode%d',tag,ib), ...
            max(abs(p{ib}.pressure-run.fields{ib}.pressure),[],'all'),0);
    end
    folder=fullfile(input_dir,[tag '_mode1']);
    fid=fopen(fullfile(folder,'case.ssp'),'r'); nr=fscanf(fid,'%d',1);
    r=1000*fscanf(fid,'%f',nr); c=fscanf(fid,'%f',nr); fclose(fid);
    wall=find(r==cfg.wall_range_m,1); assert(~isempty(wall));
    cw=c(wall); omega=2*pi*cfg.frequency_hz;
    % Discrete-Q Snell integral has a cancellation-free exact linear-segment
    % formula, independent of the continuous quadrature in the run script.
    factor=zeros(size(p{1}.pressure));
    for iz=1:size(factor,1)
        q=p{1}.receiver_depth_m(iz)-cfg.transverse_origin_m;
        for ir=1:size(factor,2)
            R=cfg.receiver_range_m(ir);
            nodes=[r(r<R);R]; cc=interp1(r,c,nodes);
            h=diff(nodes); c0=cc(1:end-1); c1=cc(2:end);
            bound=sin(deg2rad(max(abs(cfg.angle_limits_deg))))/c(1);
            fun=@(s)sum(s*(c0+c1).*h./(sqrt(1-(s*c0).^2)+sqrt(1-(s*c1).^2)))-q;
            if q==0, slow=0; else, slow=fzero(fun,[-bound bound]); end
            if strcmp(tag,'surface')
                coeff=-1;
            else
                alpha=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
                kb=omega/cfg.bottom_cp_mps+1i*alpha;
                qw=sqrt((omega/cw)^2-(omega*slow)^2);
                qb=sqrt(kb^2-(omega*slow)^2);
                if imag(qb)<0, qb=-qb; end
                coeff=(qw-qb/cfg.bottom_rho_gcm3)/(qw+qb/cfg.bottom_rho_gcm3);
            end
            factor(iz,ir)=conj(coeff);
        end
    end
    expected=p{2}.pressure.*factor;
    relative=norm(p{1}.pressure(:)-expected(:))/norm(expected(:));
    add([tag '_discrete_Q_reflected_field_relative'],relative,3e-4);
    add([tag '_exact_wall_field_nonzero'],double(any(abs(p{1}.pressure(:,1))==0)),0);
    add([tag '_root_closure_m'],max(abs(run.diag.pre_snap_wall_residual_m)),1e-8);
    add([tag '_closed_wall_residual_m'],max(abs(run.diag.wall_residual)),1e-12);
    if strcmp(tag,'surface')
        surface_residual=max(abs(p{1}.pressure(:,1)+p{2}.pressure(:,1)))/max(abs(p{2}.pressure(:,1)));
        add('surface_pressure_release_relative',surface_residual,3e-4);
    end
    reconstructed{id}=struct('cfg',cfg,'factor',factor,'field_relative',relative);
end
receipt=struct('status','COMPLETE','failed_checks',nnz(~checks.pass),'elapsed_s',toc(t), ...
    'scope','Raw SHD/MAT, exact linear-Q Snell reflection consistency, root closure; not independent FFP or total-field qualification', ...
    'matlab_release',version,'result_dir',output_dir);
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'checks','receipt','bindings','reconstructed');
fid=fopen(fullfile(output_dir,'input_sha256.json'),'w'); fprintf(fid,'%s',jsonencode(bindings,PrettyPrint=true)); fclose(fid);
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end