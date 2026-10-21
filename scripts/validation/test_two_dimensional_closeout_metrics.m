function checks = test_two_dimensional_closeout_metrics()
%TEST_TWO_DIMENSIONAL_CLOSEOUT_METRICS Analytic arrays, no acoustic solver.
names=["identical";"amplitude_only";"phase_only";"zero_candidate"; ...
    "zero_reference";"both_zero";"near_zero_reference";"missing_sample";"identity_closure"];
pass=false(numel(names),1);
m=two_dimensional_closeout_metrics([1+1i;2-1i],[1+1i;2-1i],[1;1],[true;true]);
pass(1)=m.E==0 && all(m.abs_diff==0) && all(m.phase_rad==0);
m=two_dimensional_closeout_metrics(2+0i,1+0i,1,true);
pass(2)=abs(m.E-1)<1e-14 && m.phase_rad==0 && abs(m.mag_sq-1)<1e-14 && m.phase_sq==0;
m=two_dimensional_closeout_metrics(exp(1i*pi/3),1,1,true);
pass(3)=abs(m.mag_sq)<1e-14 && abs(m.phase_sq-1)<1e-14 && abs(m.phase_rad-pi/3)<1e-14;
m=two_dimensional_closeout_metrics(0,2,1,true);
pass(4)=m.abs_diff==2 && m.rel_diff==1 && isinf(m.mag_db) && m.mag_db<0 ...
    && isnan(m.phase_rad) && m.mag_sq==4 && m.phase_sq==0;
m=two_dimensional_closeout_metrics(2,0,1,true);
pass(5)=m.abs_diff==2 && isnan(m.rel_diff) && isnan(m.mag_db) ...
    && isnan(m.phase_rad) && isnan(m.E) && m.mag_sq==4 && m.phase_sq==0;
m=two_dimensional_closeout_metrics(0,0,1,true);
pass(6)=m.abs_diff==0 && isnan(m.rel_diff) && isnan(m.phase_rad) ...
    && isnan(m.E) && m.mag_sq==0 && m.phase_sq==0;
m=two_dimensional_closeout_metrics(2e-15,1e-15,1,true);
pass(7)=m.rel_diff==1 && isfinite(m.mag_db) && m.phase_rad==0;
m=two_dimensional_closeout_metrics([1;NaN],[1;2],[1;1],[true;true]);
pass(8)=m.n_requested==2 && m.n_valid==1 && isnan(m.abs_diff(2)) && m.E==0;
A=[1+.3i;-2+1i;.2-1.3i;0]; B=[-.3+2i;3-2i;.4+.7i;2];
m=two_dimensional_closeout_metrics(A,B,ones(4,1),true(4,1));
scale=max(1,max(m.abs_diff.^2));
pass(9)=max(abs(m.closure))<1e-13*scale && abs(m.mag_fraction+m.phase_fraction-1)<1e-13;
checks=table(names,pass,'VariableNames',{'control','pass'});
assert(all(pass),'A two-dimensional closeout metric control failed.');
end
