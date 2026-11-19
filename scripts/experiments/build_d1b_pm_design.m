function receipt = build_d1b_pm_design(package_root, output_dir)
%BUILD_D1B_PM_DESIGN Surface-only frozen 2-by-2 design and direct receiver grid.
arguments
    package_root (1,1) string
    output_dir (1,1) string
end
t0 = tic;
assert(~isfolder(output_dir), 'Use a new PM audit directory');
mkdir(output_dir);
legacy = readtable(fullfile(package_root,'inputs','k0_n40.csv'));
anchor = readtable(fullfile(package_root,'inputs','k0_n12.csv'));
assert(isequal(anchor,legacy(1:12,:)) && height(legacy)==40);
L = 640; dk = 2*pi/L; m = (1:160)'; k = m*dk;
A = .74*9.81^2/6^4; Sk = .0081./(2*k.^3).*exp(-A./k.^2);
shared = 4:4:160; intermediate = find(mod(m,4)~=0 & m>=4); low = (1:3)';
assert(max(abs(k(shared)-legacy.k_rad_per_m))<1e-14);
assert(max(abs(Sk(shared)./legacy.Sk_m3-1))<1e-11);
original = [legacy.gaussian_cos legacy.gaussian_sin];
original_amplitude = [legacy.cos_depth_coeff_m legacy.sin_depth_coeff_m];
assert(max(abs(original_amplitude./sqrt(legacy.Sk_m3*4*dk)-original),[],'all')<1e-11);
Z = zeros(160,2); Z(shared,:) = original;
mid_stream = RandStream('mt19937ar','Seed',2260001);
low_stream = RandStream('mt19937ar','Seed',3260001);
Z(intermediate,:) = randn(mid_stream,numel(intermediate),2);
Z(low,:) = randn(low_stream,numel(low),2);
amp = sqrt(Sk*dk).*Z;
assert(max(abs(amp(shared,:)-original_amplitude/2),[],'all')<1e-14);
line_x = (-1139:1139)'*(200/1139); line_x([1 end]) = [-200;200];
line_z = 6*ones(size(line_x));
assert(isequal(line_x,-flipud(line_x)) && numel(unique(line_x))==2279);
write_d1_csv(table(line_x,line_z,'VariableNames',{'x_m','z_m'}),char(fullfile(output_dir,'receivers.csv')));
geometry = local_geometry(line_x,line_z);
write_d1_csv(geometry,char(fullfile(output_dir,'flat_geometry.csv')));
cases = [48 0;160 0;48 1;160 1];
rows = cell(4,1);
for j = 1:size(cases,1)
    cutoff = cases(j,1); include_low = logical(cases(j,2));
    enabled = m<=cutoff & (m>=4 | include_low);
    profile = struct('k',k(enabled),'a',amp(enabled,1),'b',amp(enabled,2),'plateau',220,'support',240);
    id = sprintf('N%03d_low%d',cutoff,include_low);
    folder = fullfile(output_dir,id); mkdir(folder);
    source_band = double(mod(m,4)==0)+2*double(m<4); % 0 intermediate, 1 shared, 2 low
    write_d1_csv(table(m,k,Sk,Z(:,1),Z(:,2),amp(:,1),amp(:,2), ...
        double(enabled),source_band,'VariableNames',{'mode','k_rad_per_m','Sk_m3','gaussian_cos', ...
        'gaussian_sin','cos_depth_m','sin_depth_m','enabled','source_band'}),char(fullfile(folder,'coefficients.csv')));
    rows{j} = local_audit(profile,folder,L,dk,Sk(enabled),line_x);
    rows{j}.id = id;
    write_d1_csv(local_stationary(profile),char(fullfile(folder,'rough_stationary_diagnostic.csv')));
    save(fullfile(folder,'profile.mat'),'profile','line_x','line_z','-v7');
end
receipt = struct('status','PASS_SURFACE_ONLY','failed_checks',0,'elapsed_s',toc(t0), ...
    'U_PM_mps',6,'seed',260001,'intermediate_seed',2260001,'low_seed',3260001, ...
    'period_m',L,'receiver_count',numel(line_x),'receiver_dx_m',200/1139, ...
    'minimum_flat_fresnel_to_platform_m',min(220-max(abs([geometry.fresnel_left_m geometry.fresnel_right_m]),[],2)), ...
    'audit',[rows{:}],'formal_cutoff_selected',false,'sound_field_solves',0);
local_json(fullfile(output_dir,'surface_receipt.json'),receipt);
end

function t = local_stationary(p)
% Sample sign changes of the geometric path derivative. Tangencies can be missed.
grid = (-9600:9600)'*.025;
[z,slope] = d1b_pm_surface(p,grid);
data = zeros(0,5);
for rx = [-200 -100 0 100 200]
    v = local_path_gradient(grid,z,slope,rx);
    brackets = find(v(1:end-1).*v(2:end)<0);
    roots_x = grid(v==0);
    for i = brackets.'
        roots_x(end+1,1) = fzero(@(x) local_gradient(p,x,rx),grid([i i+1])); %#ok<AGROW>
    end
    s = unique(roots_x); height = d1b_pm_surface(p,s);
    excess = hypot(s,1490-height)+hypot(s-rx,6-height)-hypot(rx,1496);
    data = [data;rx+zeros(size(s)) s height excess 220-abs(s)]; %#ok<AGROW>
end
t = array2table(data,'VariableNames',{'receiver_x_m','rough_stationary_x_m','wall_z_m', ...
    'path_excess_vs_flat_m','distance_to_platform_end_m'});
end

function v = local_gradient(p,x,rx)
[z,zp] = d1b_pm_surface(p,x);
v = local_path_gradient(x,z,zp,rx);
end

function v = local_path_gradient(x,z,zp,rx)
v = (x-(1490-z).*zp)./hypot(x,1490-z) + (x-rx-(6-z).*zp)./hypot(x-rx,6-z);
end

function t = local_geometry(x,z)
spec = 1490*x./(1490+z); optimal = hypot(x,1490+z);
left = zeros(size(x)); right = left;
for j = 1:numel(x)
    delta = @(s) hypot(s,1490)+hypot(x(j)-s,z(j))-optimal(j)-.375/2;
    left(j) = fzero(delta,[spec(j)-8 spec(j)]);
    right(j) = fzero(delta,[spec(j) spec(j)+8]);
end
t = table(x,z,spec,left,right,'VariableNames',{'x_m','z_m','flat_specular_x_m', ...
    'fresnel_left_m','fresnel_right_m'});
end

function metrics = local_audit(p,folder,L,dk,S,line_x)
N = 8192; x = (-N/2:N/2-1)'*(L/N);
[z,zp,zpp,raw,rp,rpp] = d1b_pm_surface(p,x);
y = [-raw -rp -rpp -z -zp -zpp]; % elevation positive up
q = fft(y)/N;
P = abs(q(1:N/2+1,:)).^2;
P(2:end-1,:) = 2*P(2:end-1,:);
moments = mean(y.^2,1);
parseval = max(abs(sum(P,1)-moments)./moments);
powers = .5*(p.a.^2+p.b.^2);
expected = sum(powers.*[ones(size(p.k)) p.k.^2 p.k.^4],1);
assert(parseval<1e-12 && max(abs(moments(1:3)./expected-1))<1e-12);
bins = round(p.k/dk)+1;
raw_leak = sum(P(setdiff(1:N/2+1,bins),1));
assert(max(abs(P(bins,1)-powers))<1e-13 && raw_leak<1e-22);
xx = (-N:N-1)'*(L/(2*N));
[zf,zpf,zppf] = d1b_pm_surface(p,xx);
fine = mean([zf zpf zppf].^2,1);
sampling = max(abs(fine./moments(4:6)-1));
assert(sampling<1e-6);
[b,bp,bpp,r,rp0,rpp0] = d1b_pm_surface(p,[-240;-220;220;240]);
assert(isequal([b bp bpp],[0 0 0;r(2:3) rp0(2:3) rpp0(2:3);0 0 0]));
outer = abs(x)>=p.support; flat = [z(outer) zp(outer) zpp(outer)];
assert(all(flat==0,'all') && all(isfinite(y),'all'));
edge = abs(x)>p.plateau & abs(x)<p.support;
% A derivative bound covers the space between samples for the clearance test.
lip = sum(p.k.*hypot(p.a,p.b))+1.875/(p.support-p.plateau)*sum(hypot(p.a,p.b));
min_wall = min(zf)-lip*L/(4*N); max_wall = max(zf)+lip*L/(4*N);
assert(min_wall>-2 && max_wall<6 && all(d1b_pm_surface(p,line_x)<6));
curvature = -zpp./(1+zp.^2).^1.5;
k = (0:N/2)'*dk;
write_d1_csv(table(k,P(:,1),P(:,4),P(:,1)/dk,P(:,4)/dk, ...
    'VariableNames',{'k_rad_per_m','raw_bin_m2','C2_bin_m2','raw_PSD_m3','C2_PSD_m3'}),char(fullfile(folder,'spectra.csv')));
write_d1_csv(array2table([x y curvature],'VariableNames', ...
    {'x_m','raw_eta_m','raw_slope','raw_second_per_m','C2_eta_m','C2_slope','C2_second_per_m','C2_curvature_per_m'}),char(fullfile(folder,'samples.csv')));
metrics = struct('enabled_modes',numel(p.k),'Kmax_rad_per_m',max(p.k), ...
    'theoretical_discrete_rms',[sqrt(dk*sum(S)) sqrt(dk*sum(S.*p.k.^2)) sqrt(dk*sum(S.*p.k.^4))], ...
    'raw_rms',sqrt(moments(1:3)),'C2_rms',sqrt(moments(4:6)), ...
    'raw_max_abs',max(abs(y(:,1:3)),[],1),'C2_max_abs',max(abs(y(:,4:6)),[],1), ...
    'C2_rms_geometric_curvature_per_m',sqrt(mean(curvature.^2)), ...
    'left_taper_rms',sqrt(mean(y(edge & x<0,4:6).^2,1)), ...
    'right_taper_rms',sqrt(mean(y(edge & x>0,4:6).^2,1)), ...
    'wall_depth_bounds_m',[min_wall max_wall],'C2_power_above_Kmax_m2',sum(P(k>max(p.k)+dk/2,4)), ...
    'raw_out_of_mode_power_m2',raw_leak,'relative_parseval_error',parseval,'relative_sample_drift',sampling);
end

function local_json(path,value)
fid = fopen(path,'w'); assert(fid>0); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
