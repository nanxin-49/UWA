function receipt = audit_bellhop_rotated_multifield(input_dirs,output_dir,options)
%AUDIT_BELLHOP_ROTATED_MULTIFIELD Independent raw per-wall and exact-Q audit.
arguments
    input_dirs (1,:) cell
    output_dir (1,:) char
    options.case_indices (1,:) cell = {}
end
assert(~isfolder(output_dir)); mkdir(output_dir); started=tic;
diary(fullfile(output_dir,'diary.log')); cleaner=onCleanup(@()diary('off'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
bindings=struct('path',{},'sha256',{}); details=cell(0,1);
for k=1:numel(input_dirs)
    data=load(fullfile(input_dirs{k},'result.mat'),'runs','receipt');
    assert(strcmp(data.receipt.status,'COMPLETE'));
    indices=1:numel(data.runs);
    if ~isempty(options.case_indices), indices=options.case_indices{k}; end
    for j=indices
        if isempty(data.runs{j}), continue; end
        a=data.runs{j}; cfg=a.cfg; N=cfg.reflection_order;
        tag=sprintf('%s_N%d',cfg.direction,N);
        folder=fullfile(input_dirs{k},sprintf('%s_order%d_mode1',cfg.direction,N));
        stem=fullfile(folder,'case'); raw=read_bellhop_rotated_field(stem);
        for ext={'.env','.ssp','.iwcov','.rotq','.r25d','.shd','.prt','.covdiag','.mwdiag','.r25diag'}
            path=[stem ext{1}]; bindings(end+1)=struct('path',path,'sha256',black_file_hash(path)); %#ok<AGROW>
        end
        add([tag '_raw_SHD_binding'],max(abs(raw.pressure-a.fields{1}.pressure),[],'all'),0);
        physical=complex(zeros(size(a.physical_pressure)));
        physical(cfg.receiver_depth_order,cfg.receiver_x_order)=(-conj(raw.pressure)).';
        add([tag '_physical_permutation'],max(abs(physical-a.physical_pressure),[],'all'),0);
        fid=fopen([stem '.ssp'],'r'); nr=fscanf(fid,'%d',1); r=1000*fscanf(fid,'%f',nr);
        c=fscanf(fid,'%f',nr); c2=fscanf(fid,'%f',nr); fclose(fid);
        add([tag '_transverse_constant_Q'],max(abs(c2-c)),0);
        sg=1-2*strcmp(cfg.direction,'up'); z0=cfg.source_z_m+sg*r(1:end-1);
        z=cfg.water_depth_m-abs(mod(z0,2*cfg.water_depth_m)-cfg.water_depth_m);
        add([tag '_folded_Munk_Q'],max(abs(c(1:end-1)-munk_profile_sound_speed(z,cfg.munk_parameters))),2e-10);
        wall=readmatrix([stem '.covdiag'],'FileType','text','CommentStyle','#'); wall=wall(~all(isnan(wall),2),:);
        assert(size(wall,2)==44);
        add([tag '_all_wall_records'],abs(height(wall)-N*cfg.beam_count),0);
        wall_number=mod((0:height(wall)-1).',N)+1;
        expected_r=cfg.wall_range_m+(wall_number-1)*cfg.water_depth_m;
        add([tag '_all_wall_locations_m'],max(abs(wall(:,3)-expected_r)),1e-10);
        add([tag '_closed_wall_residual_m'],max(abs(wall(:,5))),1e-12);
        add([tag '_proper_rotation'],max(wall(:,17)),1e-12);
        add([tag '_specular_direction'],max(wall(:,16)),1e-8);
        add([tag '_q_preserved_at_reflection'],max(abs(wall(:,33))),1e-10);
        add([tag '_pq_preserved_by_rotation'],max(abs(wall(:,34:35)),[],'all'),1e-10);
        surface=(mod(wall_number,2)==1)==strcmp(cfg.direction,'up');
        omega=2*pi*cfg.frequency_hz;
        cw=munk_profile_sound_speed(cfg.water_depth_m,cfg.munk_parameters);
        alpha=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
        kb=omega/cfg.bottom_cp_mps+1i*alpha;
        sin_inc=wall(:,11)./hypot(wall(:,10),wall(:,11));
        kr=(omega/cw)*sin_inc;
        qw=sqrt((omega/cw)^2-kr.^2); qb=sqrt(kb^2-kr.^2);
        qb(imag(qb)<0)=-qb(imag(qb)<0);
        coeff=conj((qw-qb/cfg.bottom_rho_gcm3)./(qw+qb/cfg.bottom_rho_gcm3)); coeff(surface)=-1;
        observed=(wall(:,22)./wall(:,21)).*exp(1i*(wall(:,19)-wall(:,18)));
        add([tag '_each_complex_reflection'],max(abs(observed-coeff)),1e-10);
        d=readmatrix([stem '.r25diag'],'FileType','text','CommentStyle','#'); d=d(~all(isnan(d),2),:);
        time=zeros(height(d),1); integ=time;
        for b=1:height(d)
            R=d(b,2); rn=[r(r<R);R]; cn=interp1(r,c,rn);
            h=diff(rn); ca=cn(1:end-1); cb=cn(2:end); dc=cb-ca;
            p=sind(d(b,1))/c(1); sa=sqrt(1-(p*ca).^2); sb=sqrt(1-(p*cb).^2);
            dt=h./(ca.*sa); varying=abs(dc)>1e-12;
            dt(varying)=h(varying)./dc(varying).*(log1p(dc(varying)./ca(varying))+ ...
                log1p((sa(varying)-sb(varying))./(1+sb(varying))));
            time(b)=sum(dt); integ(b)=sum(h.*(ca+cb)./(sa+sb));
        end
        add([tag '_exact_linear_Q_phase_rad'],omega*max(abs(d(:,4)-time)),.01);
        add([tag '_exact_outplane_integral_relative'],max(abs(d(:,6)-integ)./integ),1e-7);
        details{end+1}=struct('tag',tag,'exact_time_s',time,'exact_integral_cds',integ); %#ok<AGROW>
    end
end
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status','COMPLETE','failed_checks',nnz(~checks.pass),'elapsed_s',toc(started), ...
    'scope','Raw per-wall complex coefficients, reflection geometry, exact linear-Q time/integral and physical permutation; no complete-field or FFP qualification', ...
    'input_dirs',{input_dirs},'case_indices',{options.case_indices},'result_dir',output_dir);
save(fullfile(output_dir,'result.mat'),'receipt','checks','details','bindings');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
write(fullfile(output_dir,'input_sha256.json'),jsonencode(bindings,PrettyPrint=true));
disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end