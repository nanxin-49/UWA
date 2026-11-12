function receipt=validate_munk_ffp_long_uniform(run_id)
%VALIDATE_MUNK_FFP_LONG_UNIFORM Uniform limit at the actual 5000m/4950m scale.
% Existing point/image references reused; full spectral kernel integrated with
% analytic subtraction disabled, and error checked at EACH receiver depth.
arguments
    run_id (1,:) char = '20261112_full_01'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'scripts','validation','support'),'-begin');
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id);
assert(~isfile(fullfile(d,'long_uniform_receipt.json'))); t=tic; diary(fullfile(d,'long_uniform_diary.log'));
receipt=struct('stage','long_uniform_limit','status','RUNNING','elapsed_s',0,'exit_status',NaN);
try
    cfg=make_munk_5000m_ffp_config(5000); cfg.layer_c_mps=1500;
    cfg.depth_m=[0;1000;3000;4900;4948;4952;4999;5000]; cfg.range_m=[-5 -.1 0 .1 5];
    cfg.bottom_type='outgoing'; cfg.top_type='outgoing'; cfg.subtract_reference=false;
    cfg.panel_width_radpm=.0625; cfg.kappa_max_radpm=100; cfg.contour_offset_radpm=.1; cfg.spectral_block=256; cfg.max_wall_s=90;
    free=solve_layered_fluid_ffp(cfg);
    cfg.top_type='pressure_release'; release=solve_layered_fluid_ffp(cfg);
    [zz,xx]=ndgrid(cfg.depth_m,cfg.range_m); pts=[xx(:),zeros(numel(xx),1),zz(:)];
    exact_free=reshape(deepsea3d_point_field(pts,[0 0 4950],2*pi*1e4/1500),size(xx));
    exact_release=complex(zeros(size(exact_free)));
    for j=1:numel(exact_release)
        if pts(j,3)==0, continue; end % exact Dirichlet boundary, existing image helper requires z>0
        p=struct('frequency_hz',1e4,'c0_mps',1500,'tx_m',[0 0 4950],'rx_m',pts(j,:));
        old=deepsea3d_flat(p,0); exact_release(j)=old.total_native;
    end
    rows=struct([]);
    for m=1:2
        if m==1, actual=free.green; expected=exact_free; name='free'; else, actual=release.green; expected=exact_release; name='surface_mirror'; end
        for j=1:numel(cfg.depth_m)
            den=norm(expected(j,:));
            if den==0
                rel=NaN; tl=NaN; phase=NaN; pass=max(abs(actual(j,:)))<=1e-12;
            else
                rel=norm(actual(j,:)-expected(j,:))/den;
                tl=max(abs(20*log10(abs(actual(j,:))./abs(expected(j,:)))));
                phase=max(abs(angle(actual(j,:).*conj(expected(j,:)))));
                pass=rel<=1e-3 && tl<=.05 && phase<=.01;
            end
            rows=[rows struct('reference',string(name),'depth_m',cfg.depth_m(j),'complex_rel_l2',rel, ...
                'max_tl_error_db',tl,'max_phase_error_rad',phase,'max_abs_error',max(abs(actual(j,:)-expected(j,:))),'pass',pass)]; %#ok<AGROW>
        end
    end
    save(fullfile(d,'long_uniform_controls.mat'),'free','release','exact_free','exact_release','rows','-v7');
    tbl=struct2table(rows); writetable(tbl,fullfile(d,'long_uniform_checks.csv'));
    receipt.failed_checks=nnz(~tbl.pass); receipt.all_depths_pass=all(tbl.pass);
    finite=isfinite(tbl.complex_rel_l2);
    receipt.metrics=struct('max_per_depth_complex_rel',max(tbl.complex_rel_l2(finite)), ...
        'max_tl_error_db',max(tbl.max_tl_error_db(finite)),'max_phase_error_rad',max(tbl.max_phase_error_rad(finite)), ...
        'receiver_depths_m',cfg.depth_m.','source_depth_m',4950,'analytic_subtraction_used',false);
    if receipt.all_depths_pass, receipt.status='ACTUAL_DEPTH_SCALE_UNIFORM_LIMIT_PASS'; receipt.exit_status=0;
    else, receipt.status='ACTUAL_DEPTH_SCALE_UNIFORM_LIMIT_FAIL'; receipt.exit_status=1; end
catch err
    receipt.status='FAILED'; receipt.exit_status=1; receipt.error=getReport(err,'extended','hyperlinks','off');
end
receipt.elapsed_s=toc(t); fid=fopen(fullfile(d,'long_uniform_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,'PrettyPrint',true)); fclose(fid); diary off;
fprintf('%s\n',jsonencode(receipt));
end
