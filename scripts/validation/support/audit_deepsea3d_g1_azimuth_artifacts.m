function checks=audit_deepsea3d_g1_azimuth_artifacts()
%AUDIT_DEEPSEA3D_G1_AZIMUTH_ARTIFACTS Verify fixed sea, native fields and contrast.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_g1_azimuth');
a=load(fullfile(out,'azimuth.mat'),'state');s=a.state;
assert(s.complete && s.global_G1=="OPEN" && ~s.WP2_allowed);
assert(isequal(s.cfg.tx_m,[0 150.5 1500]) && isequal(s.baseline.cfg.tx_m,[150.5 0 1500]));
assert(isequal(s.cfg.rx_m,[0 0 5]) && norm(s.cfg.tx_m(1:2))==norm(s.baseline.cfg.tx_m(1:2)));
restored=s.cfg;restored.tx_m=s.baseline.cfg.tx_m;assert(isequal(restored,s.baseline.cfg));
assert(s.flat_control.passed && s.flat_control.x_y_direct_error<1e-10 && s.flat_control.x_y_reflect_error<1e-10);
x=-12:.0625:12;[xx,yy]=meshgrid(x,x);
[e0,dx0,dy0]=deepsea3d_surface(xx,yy,s.baseline.cfg.profile);
[e1,dx1,dy1]=deepsea3d_surface(xx,yy,s.cfg.profile);
surface_error=max(abs([e0(:)-e1(:);dx0(:)-dx1(:);dy0(:)-dy1(:)]));assert(surface_error==0);
checks=struct('complete',true,'azimuth_only_change',true,'fixed_surface_max_difference',surface_error, ...
    'flat_backgrounds_equal',true,'reference_ready',s.reference_ready, ...
    'azimuth_G1',s.G1_status,'global_G1',"OPEN",'WP2_allowed',false);
if ~s.reference_ready
    assert(isempty(s.fk) && s.G1_status=="INCONCLUSIVE");
else
    g=readtable(fullfile(out,'reference_guards.csv'));assert(all(g{1,:}<1e-4));
    r=s.runs(s.reference_index).result;o=s.observer;
    assert(all(isfinite([r.terms;r.partials;r.boundary_coefficient_norms])));
    assert(all(diff(r.boundary_coefficient_norms(end-4:end))<0) && all(diff(abs(r.terms(end-4:end)))<0));
    assert(o.term_match_error<1e-12 && o.norm_match_error<1e-12 && r.forcing_series_relative<1e-8);
    assert(norm(o.terms-r.terms)/norm(r.terms)<1e-12);
    f=readtable(fullfile(out,'fk_runs.csv'));refs=readtable(fullfile(out,'reference_runs.csv'));
    assert(height(f)==5 && height(refs)==numel(s.runs));
    F=f.reflect_real(2)+1i*f.reflect_imag(2);D=f.delta_real(2)+1i*f.delta_imag(2);
    R=refs.reflect_real(s.reference_index)+1i*refs.reflect_imag(s.reference_index);
    I=refs.delta_real(s.reference_index)+1i*refs.delta_imag(s.reference_index);
    assert(abs(F-s.fk(2).reflect_native)/abs(F)<1e-13 && abs(I-r.rough_correction)/abs(I)<1e-13);
    ratio=[F/R,D/I];errors=abs(ratio-1);phase=angle(ratio);rho=abs(ratio);
    identity=max(abs(errors.^2-((rho-1).^2+4*rho.*sin(phase/2).^2)));assert(identity<1e-16);
    c=summarize_deepsea3d_azimuth_comparison(s);assert(isequaln(c,s.comparison));
    p=readtable(fullfile(out,'azimuth_comparison.csv'));
    assert(max(abs(errors-[p.reflected_complex(2),p.increment_complex(2)]))<1e-13);
    assert(max(abs(phase-[p.reflected_phase_rad(2),p.increment_phase_rad(2)]))<1e-13);
    old=readtable(fullfile(root,'results','validation','deepsea3d_g1_q10','geometry_comparison.csv'));
    assert(abs(p.reflected_complex(1)-old.reflected_complex(2))<1e-14);
    assert(abs(p.increment_complex(1)-old.increment_complex(2))<1e-14);
    contrast=readtable(fullfile(out,'directionality.csv'));
    native_contrast=R-complex(p.reference_reflect_real(1),p.reference_reflect_imag(1));
    assert(abs(native_contrast-complex(contrast.y_minus_x_real(1),contrast.y_minus_x_imag(1)))<1e-18);
    assert(abs((R-I)-(complex(p.reference_reflect_real(1),p.reference_reflect_imag(1))- ...
        complex(p.reference_delta_real(1),p.reference_delta_imag(1))))<1e-18);
    num=s.fk_guards;
    assert(classify_deepsea3d_amplitude(true,max([num.grid,num.aperture_12_16,num.aperture_16_24]),num.gauss,c.metrics_y)==s.G1_status);
    sensitivity=local_order_sensitivity(s,out);
    checks.CSV_native_pressure_phase_matches=true;checks.amplitude_phase_identity_abs_error=identity;
    checks.x_case_saved_metrics_unchanged=true;checks.reference_directionality_recomputed=true;
    checks.order_sensitivity_variants_each=sensitivity.reference_variants_each;
end
local_json(fullfile(out,'saved_data_checks.json'),checks);disp(checks);
end
function diagnostic=local_order_sensitivity(s,out)
% Post-result supplement: exact lower-order pressure partials, no new solve.
% Lower-order boundary residuals are not separately certified by this replay.
rx=local_pressure_variants(s.baseline);ry=local_pressure_variants(s);
fields={'reflect_native','rough_correction'};rows=struct([]);
for j=1:2
    field=fields{j};vx=local_errors(s.baseline.fk,rx,field);vy=local_errors(s.fk,ry,field);
    row=struct('field',string(field),'x_min',min(vx),'x_max',max(vx),'y_min',min(vy),'y_max',max(vy), ...
        'ranges_disjoint',max(vy)<min(vx) || min(vy)>max(vx));
    if j==1,rows=repmat(row,2,1);else,rows(j)=row;end
end
main_x=s.baseline.runs(s.baseline.reference_index).result;main_y=s.runs(s.reference_index).result;
contrast=main_y.rough_correction-main_x.rough_correction;spread=0;
for x=rx
    for y=ry
        spread=max(spread,abs((y.rough_correction-x.rough_correction)-contrast));
    end
end
diagnostic=struct('reference_variants_each',numel(rx),'FK_variants_each',5,'envelopes',rows, ...
    'contrast_abs',abs(contrast),'contrast_max_abs_displacement',spread, ...
    'contrast_spread_relative',spread/abs(contrast),'contrast_resolved',abs(contrast)>spread, ...
    'post_result_supplement',true,'scope','selected N, N-2, N-4 pressure partials on three grids; lower-order boundary fields not independently certified; sensitivity only');
writetable(struct2table(rows),fullfile(out,'order_sensitivity_envelopes.csv'));
local_json(fullfile(out,'order_sensitivity.json'),diagnostic);disp(diagnostic);
end
function variants=local_pressure_variants(s)
indices=[s.coarse_index,s.fine_index,s.reference_index];variants=struct([]);bg=deepsea3d_flat(s.cfg,0);
for index=indices
    r=s.runs(index).result;
    for N=[r.order-4,r.order-2,r.order]
        d=r.partials(N);p=struct('reflect_native',bg.reflect_native+d,'rough_correction',d);
        if isempty(variants),variants=p;else,variants(end+1)=p;end %#ok<AGROW>
    end
end
end
function values=local_errors(fk,refs,field)
values=zeros(numel(fk)*numel(refs),1);i=0;
for f=fk
    for r=refs
        i=i+1;values(i)=abs(f.(field)-r.(field))/abs(r.(field));
    end
end
end
function local_json(file,s)
fid=fopen(file,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(s,PrettyPrint=true));
end
