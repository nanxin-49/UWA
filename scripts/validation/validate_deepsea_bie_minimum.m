function output=validate_deepsea_bie_minimum()
%VALIDATE_DEEPSEA_BIE_MINIMUM Matrix-free reference feasibility, no PM solve.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','deepsea_bie_reconstruction_20261014');
file=fullfile(outdir,'bie_minimum.mat');assert(~isfile(file),'BIE:Exists','No overwrite.');
totalclock=tic;[before,sys]=memory;reserve=.3*sys.PhysicalMemory.Available/2^20;
assert(384<.7*sys.PhysicalMemory.Available/2^20,'BIE:RAM','384MiB estimate fails physical reserve.');
proc=System.Diagnostics.Process.GetCurrentProcess;proc.Refresh();
peak0=double(proc.PeakWorkingSet64)/2^20;working0=double(proc.WorkingSet64)/2^20;
samples=[0,sys.PhysicalMemory.Available/2^20,working0];lastsample=tic;
tests=struct();
[xy,w]=deepsea_bie_duffy([0 0],1,[0 0],32);
points=[xy,zeros(size(xy,1),1)];normal=repmat([0 0 1],size(xy,1),1);
[S,~]=deepsea_bie_layers([0 0 0],points,normal,w,ones(size(w)),0,Inf);
tests.Laplace_self_relative=abs(S-asinh(1)/pi)/(asinh(1)/pi);
assert(tests.Laplace_self_relative<1e-10);
near=struct([]);
for degree=[8 16 24]
    [xy,w]=deepsea_bie_duffy([0 0],1,[0 0],degree);ps=[xy,zeros(size(xy,1),1)];
    ns=repmat([0 0 1],numel(w),1);
    for epsilon=[.1 .01 .001 -.001]
        [~,Dk]=deepsea_bie_layers([0 0 epsilon],ps,ns,w,ones(size(w)),2.5,Inf);
        [~,D0]=deepsea_bie_layers([0 0 epsilon],ps,ns,w,ones(size(w)),0,Inf);
        solid=sign(epsilon)/pi*atan(.25/(abs(epsilon)*sqrt(epsilon^2+.5)));
        corrected=Dk-D0+solid;
        row=struct('degree',degree,'epsilon',epsilon,'raw_D_real',real(Dk), ...
            'subtracted_D_real',real(corrected),'subtracted_D_imag',imag(corrected), ...
            'solid_angle',solid,'jump_limit_error',abs(corrected-sign(epsilon)*.5));
        if isempty(near),near=row;else,near(end+1)=row;end %#ok<AGROW>
    end
end
tests.jump_sign=near(end).subtracted_D_real<0 && near(end-1).subtracted_D_real>0;
assert(tests.jump_sign);
% Image-normal derivative checked without using the implemented D kernel.
target=[.4 .1 .8];source=[.1 -.3 .03];nn=[.1 .2 1];nn=nn/norm(nn);e=1e-5;
[~,D]=deepsea_bie_layers(target,source,nn,1,1,2.5,-1);
v=zeros(1,2);
for j=1:2
    y=source+(2*j-3)*e*nn;im=y;im(3)=-2-y(3);
    R=norm(target-y);Ri=norm(target-im);v(j)=exp(2.5i*R)/(4*pi*R)-exp(2.5i*Ri)/(4*pi*Ri);
end
tests.image_normal_FD_relative=abs((v(2)-v(1))/(2*e)-D)/abs(D);
assert(tests.image_normal_FD_relative<1e-7);
raw={};cases=struct([]);caseclock=tic;calls=0;physical_kernel=0;
for frequency=[600 6000]
    for n=[12 20 40]
        caseclock=tic;calls=0;guard();k=2*pi*frequency/1500;
        p=deepsea_bie_flat_plan(3,n,k,-1);sigma=.7;
        rhs=deepsea_bie_gaussian_reference(p.target,k,-1,sigma,256);
        rhs_fine=deepsea_bie_gaussian_reference(p.target,k,-1,sigma,512);
        reference_change=norm(rhs-rhs_fine)/norm(rhs_fine);
        % Right-hand side uses the independent finer Fourier reference.
        try
            [phi,flag,relres,iter,hist]=gmres(@action,rhs_fine,8,1e-10,4);
            fresh=norm(action(phi)-rhs_fine)/norm(rhs_fine);
            shifted=deepsea_bie_flat_plan(3,n,k,-1,[.27 .19]);
            boundary=deepsea_bie_flat_action(phi,shifted,@guard);
            truth_boundary=deepsea_bie_gaussian_reference(shifted.target,k,-1,sigma,512);
            boundary_error=norm(boundary-truth_boundary)/norm(truth_boundary);
            obs=[0 0 .3;.5 .2 1;0 0 5];
            [Sf,Df]=deepsea_bie_layers(obs,p.source,p.normal,p.weight,phi(p.cell_index),k,-1,@guard);
            field=Df-1i*k*Sf;truth_field=deepsea_bie_gaussian_reference(obs,k,-1,sigma,512);
            field_error=norm(field-truth_field)/norm(truth_field);
            density_exact=exp(-sum(p.center.^2,2)/(2*sigma^2));
            row=struct('frequency_hz',frequency,'n',n,'DOF',n^2,'h_m',p.h, ...
                'gmres_flag',flag,'gmres_relative',relres,'fresh_linear',fresh, ...
                'independent_boundary',boundary_error,'independent_field',field_error, ...
                'density_error',norm(phi-density_exact)/norm(density_exact), ...
                'reference_change',reference_change,'operator_calls',calls,'inner_iterations',(iter(1)-1)*8+iter(2), ...
                'linear_pass',fresh<1e-8,'reference_pass',max([boundary_error,field_error])<1e-4 && reference_change<1e-6, ...
                'wall_s',toc(caseclock),'status',"COMPLETED");
            item=struct('phi',phi,'rhs',rhs_fine,'history',hist,'boundary',boundary, ...
                'truth_boundary',truth_boundary,'field',field,'truth_field',truth_field);
        catch err
            row=struct('frequency_hz',frequency,'n',n,'DOF',n^2,'h_m',p.h, ...
                'gmres_flag',-99,'gmres_relative',NaN,'fresh_linear',NaN, ...
                'independent_boundary',NaN,'independent_field',NaN,'density_error',NaN, ...
                'reference_change',reference_change,'operator_calls',calls,'inner_iterations',NaN, ...
                'linear_pass',false,'reference_pass',false,'wall_s',toc(caseclock),'status',string(err.identifier));
            item=struct('error',err.message);
        end
        physical_kernel=physical_kernel+row.wall_s;
        if isempty(cases),cases=row;else,cases(end+1)=row;end %#ok<AGROW>
        raw{end+1}=item; %#ok<AGROW>
        fprintf('BIE manufactured %dHz n%d linear%.3g boundary%.3g field%.3g wall%.2fs\n', ...
            frequency,n,row.fresh_linear,row.independent_boundary,row.independent_field,row.wall_s);
        if row.status~="COMPLETED",break;end
    end
end
% The known-data compact correction tests geometry, not an unknown PM solve.
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
cauchy=deepsea_bie_cauchy_control(a.state.dual_cfg.profile);
physical_kernel=physical_kernel+sum([cauchy.wall_s]);
[after,sys]=memory;proc.Refresh();
output=struct('cases',cases,'local_tests',tests,'near_jump',near,'cauchy',cauchy, ...
    'A_kernel_wall_s',physical_kernel,'total_wall_s',toc(totalclock), ...
    'memory_before_MiB',before.MemUsedMATLAB/2^20,'memory_after_MiB',after.MemUsedMATLAB/2^20, ...
    'available_after_MiB',sys.PhysicalMemory.Available/2^20,'memory_samples',samples, ...
    'lifetime_working_peak_before_MiB',peak0,'lifetime_working_peak_after_MiB',double(proc.PeakWorkingSet64)/2^20, ...
    'peak_scope','process lifetime upper bound, plus task samples; not isolated task allocation peak', ...
    'unknown_rough_exterior_closure_verified',false,'PM_BIE_solve',false, ...
    'reference_ready',false,'G1_target_PM','OPEN','WP2_formal',false);
assert(output.A_kernel_wall_s<3600 && output.total_wall_s<3600);
save(file,'output','raw','-v7');writetable(struct2table(cases),fullfile(outdir,'bie_manufactured.csv'));
writetable(struct2table(cauchy),fullfile(outdir,'bie_cauchy.csv'));writetable(struct2table(near),fullfile(outdir,'bie_near_jump.csv'));
fid=fopen(fullfile(outdir,'bie_minimum.json'),'w');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
    function z=action(v)
        calls=calls+1;if calls>80,error('BIE:Calls','80 application cap.');end
        z=deepsea_bie_flat_action(v,p,@guard);
    end
    function guard()
        if toc(totalclock)>3550 || toc(caseclock)>180,error('BIE:Time','Declared case/study cap.');end
        if toc(lastsample)>.5
            [~,s]=memory;proc.Refresh();available=s.PhysicalMemory.Available/2^20;
            samples(end+1,:)=[toc(totalclock),available,double(proc.WorkingSet64)/2^20];
            if available<reserve,error('BIE:RAM','Preserve initial 30 percent available RAM.');end
            lastsample=tic;
        end
    end
end
