function receipt = audit_piston_source_quicklook(run_dir)
%AUDIT_PISTON_SOURCE_QUICKLOOK Independent spectral/adaptive saved-data checks.
% This limited visualization audit does not implement P0 gradients/integration.
arguments
    run_dir (1,:) char
end
t0=tic; rows=[]; metrics=[]; direction_table=table();
for fk=[4 6 8]
    s=load(fullfile(run_dir,sprintf('field_%dkHz.mat',fk)));
    cfg=s.cfg; f=s.f_hz; k=2*pi*f/cfg.c_mps;
    T=readtable(fullfile(run_dir,sprintf('receiver_%dkHz.csv',fk)));
    for wi=1:3
        W=cfg.total_widths_m(wi); coeff=cfg.rho_kgpm3*2*pi*f*cfg.v0_mps/2;
        % Independent adaptive spatial integration at near/remote sentinels.
        points=[0 .01;.10 .01;-.20 .10;.13 .40;0 97;40 97;120 97];
        ph=piston_strip_pressure_quicklook(points(:,1),points(:,2),W,f,cfg,cfg.near_nquad);
        pa=complex(zeros(size(ph)));
        for jj=1:numel(ph)
            x=points(jj,1); d=points(jj,2);
            pa(jj)=coeff*integral(@(xi)besselh(0,1,k*hypot(x-xi,d)),-W/2,W/2, ...
                'RelTol',1e-11,'AbsTol',1e-12);
        end
        adaptive_rel=max(abs(ph-pa)./abs(pa));
        % Different representation: full propagating AND evanescent spectra.
        as_points=[0 .02;.13 .04;-.2 .1;.2 .4];
        ph_as=piston_strip_pressure_quicklook(as_points(:,1),as_points(:,2),W,f,cfg,cfg.near_nquad);
        ps=complex(zeros(size(ph_as)));
        for jj=1:numel(ps)
            x=as_points(jj,1); d=as_points(jj,2);
            prop=integral(@(a)local_sinc(k*W/2*sin(a)).*cos(k*x*sin(a)).*exp(1i*k*d*cos(a)), ...
                0,pi/2,'RelTol',1e-11,'AbsTol',1e-12);
            tmax=asinh(40/(k*d));
            evan=integral(@(t)local_sinc(k*W/2*cosh(t)).*cos(k*x*cosh(t)).*exp(-k*d*sinh(t)), ...
                0,tmax,'RelTol',1e-11,'AbsTol',1e-12);
            ps(jj)=cfg.rho_kgpm3*2*pi*f*cfg.v0_mps*W/pi*(prop-1i*evan);
        end
        as_rel=max(abs(ps-ph_as)./abs(ps));
        % Re-evaluate a two-dimensional lattice from both saved maps.
        [I,J]=ndgrid(1:19:numel(s.z_m),1:31:numel(s.x_m));
        idx=sub2ind(size(s.p_field_f(:,:,wi)),I(:),J(:));
        ref=piston_strip_pressure_quicklook(reshape(s.x_m(J(:)),[],1), ...
            reshape(100-s.z_m(I(:)),[],1),W,f,cfg,48);
        pp=s.p_field_f(:,:,wi); whole_rel=max(abs(pp(idx)-ref)./max(abs(ref),1e-10));
        [I,J]=ndgrid(1:17:numel(s.s_near_m),1:13:numel(s.x_near_m));
        idx=sub2ind(size(s.p_near_f(:,:,wi)),I(:),J(:));
        ref=piston_strip_pressure_quicklook(reshape(s.x_near_m(J(:)),[],1), ...
            reshape(s.s_near_m(I(:)),[],1),W,f,cfg,192);
        pp=s.p_near_f(:,:,wi); near_rel=max(abs(pp(idx)-ref)./max(abs(ref),1e-10));
        tag=sprintf('W%02dcm',round(100*W));
        pcsv=T.([tag '_real_Pa'])+1i*T.([tag '_imag_Pa']);
        csv_rel=max(abs(pcsv-s.p_rx_f(:,wi)))/max(abs(s.p_rx_f(:,wi)));
        symmetry=max(abs(s.p_rx_f(:,wi)-flipud(s.p_rx_f(:,wi))))/max(abs(s.p_rx_f(:,wi)));
        iz=find(s.z_m==3); [~,irx]=ismember(s.x_m,s.x_rx_m);
        pmap=reshape(s.p_field_f(iz,:,wi),[],1);
        line_rel=max(abs(pmap-s.p_rx_f(irx,wi)))/max(abs(pmap));
        dn=abs(s.p_circle_f(:,wi))/abs(s.p_circle_f(s.theta_deg==0,wi));
        far_abs=max(abs(dn-abs(s.D_sinc(:,wi))));
        % Pressure far-field approximation, not used for the pressure maps.
        r=hypot(s.x_rx_m,97); theta=atan2(s.x_rx_m,97);
        pff=coeff*W*sqrt(2./(pi*k*r)).*exp(1i*(k*r-pi/4)).*local_sinc(k*W/2*sin(theta));
        far_complex=max(abs(pff-s.p_rx_f(:,wi)))/max(abs(s.p_rx_f(:,wi)));
        failed=adaptive_rel>1e-7 || as_rel>1e-7 || whole_rel>1e-7 || near_rel>1e-7 || ...
            csv_rel>1e-12 || symmetry>1e-12 || line_rel>1e-12 || far_abs>1e-3 || far_complex>1e-3;
        rows=[rows;f W adaptive_rel as_rel whole_rel near_rel csv_rel symmetry line_rel far_abs far_complex failed]; %#ok<AGROW>
        pos=s.theta_deg>=0;
        tt=s.theta_deg(pos); dd=dn(pos);
        cross=find(dd<=10^(-3/20),1);
        if isempty(cross), hpbw=NaN; else
            half=interp1(dd(cross-1:cross),tt(cross-1:cross),10^(-3/20)); hpbw=2*half;
        end
        if cfg.c_mps/f/W<=1, null_deg=asind(cfg.c_mps/f/W); else, null_deg=NaN; end
        axis_pa=abs(s.p_rx_f(s.x_rx_m==0,wi));
        at40=abs(s.p_rx_f(s.x_rx_m==40,wi));
        at120=abs(s.p_rx_f(s.x_rx_m==120,wi));
        metrics=[metrics; f W W/(cfg.c_mps/f) axis_pa at40 20*log10(at40/axis_pa) at120 20*log10(at120/axis_pa) hpbw null_deg]; %#ok<AGROW>
        D=table(f+zeros(size(s.theta_deg)),W+zeros(size(s.theta_deg)),s.theta_deg, ...
            real(s.p_circle_f(:,wi)),imag(s.p_circle_f(:,wi)),dn,s.D_sinc(:,wi), ...
            'VariableNames',{'frequency_hz','total_width_m','theta_deg','real_Pa','imag_Pa','exact_axis_normalized','analytic_sinc_signed'});
        direction_table=[direction_table;D]; %#ok<AGROW>
    end
end
A=array2table(rows,'VariableNames',{'frequency_hz','total_width_m','adaptive_relative','AS_relative', ...
    'whole_map_quad_relative','near_map_quad_relative','csv_relative','symmetry_relative','map_line_relative', ...
    'farfield_normalized_absolute','farfield_complex_relative','failed'});
M=array2table(metrics,'VariableNames',{'frequency_hz','total_width_m','width_over_lambda', ...
    'axis_peak_Pa','x40_peak_Pa','x40_vs_axis_dB','x120_peak_Pa','x120_vs_axis_dB','full_minus3dB_width_deg','first_sinc_null_deg'});
writetable(A,fullfile(run_dir,'audit_metrics.csv')); writetable(M,fullfile(run_dir,'summary_metrics.csv'));
writetable(direction_table,fullfile(run_dir,'direction_curves.csv'));
receipt=struct('status','QUICKLOOK_CHECKS_COMPLETE','failed_checks',sum(A.failed), ...
    'elapsed_s',toc(t0),'max_AS_complex_relative',max(A.AS_relative), ...
    'max_map_quadrature_relative',max([A.whole_map_quad_relative;A.near_map_quad_relative]), ...
    'max_farfield_normalized_absolute',max(A.farfield_normalized_absolute), ...
    'unrun_out_of_scope',{{'P0 pressure gradients','production shared-source integration','Gaussian interface regression','P1 plane sea surface'}}, ...
    'qualification','visualization numerical checks only; no source qualification');
fid=fopen(fullfile(run_dir,'audit_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt)); assert(receipt.failed_checks==0,'Quicklook numerical audit failed; see audit_metrics.csv.');
end
function y=local_sinc(q)
y=ones(size(q)); nz=q~=0; y(nz)=sin(q(nz))./q(nz);
end
