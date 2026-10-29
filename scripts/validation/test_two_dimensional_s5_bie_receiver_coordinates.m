function tests = test_two_dimensional_s5_bie_receiver_coordinates
%TEST_TWO_DIMENSIONAL_S5_BIE_RECEIVER_COORDINATES Receiver evaluation contract.
tests = functiontests(localfunctions);
end

function testPairedCoordinatesMatchSeparateEvaluations(testCase)
cfg = local_config();
xv = [-0.1953125;0;0.1953125];
zv = [2.875;3;3.125];
cfg.receiver_x_m = repmat(xv,3,1);
cfg.receiver_z_m = repelem(zv,3);
ids = ["P11";"P12";"P13";"P21";"P22";"P23";"P31";"P32";"P33"];
for iz = 1:3
    for ix = 1:3
        idx = (iz-1)*3+ix;
        verifyEqual(testCase,ids(idx),string(sprintf('P%d%d',iz,ix)));
        verifyEqual(testCase,cfg.receiver_x_m(idx),xv(ix));
        verifyEqual(testCase,cfg.receiver_z_m(idx),zv(iz));
    end
end
paired = solve_helmholtz_bie_halfplane_vertical(cfg);

separate = complex(zeros(size(cfg.receiver_x_m)));
for ii = 1:numel(separate)
    point_cfg = cfg;
    point_cfg.receiver_x_m = cfg.receiver_x_m(ii);
    point_cfg.receiver_z_m = cfg.receiver_z_m(ii);
    point = solve_helmholtz_bie_halfplane_vertical(point_cfg);
    separate(ii) = point.receiver_field;
end

relative_error = norm(paired.receiver_field - separate) / norm(separate);
verifyLessThan(testCase, relative_error, 1e-12);
verifyEqual(testCase, paired.receiver_z_m, cfg.receiver_z_m);
verifyEqual(testCase, paired.receiver_x_m, cfg.receiver_x_m);
end

function testScalarDepthStillExpandsAcrossReceivers(testCase)
cfg = local_config();
cfg.receiver_x_m = [-0.4; 0; 0.4];
cfg.receiver_z_m = 3.0;
out = solve_helmholtz_bie_halfplane_vertical(cfg);

verifyEqual(testCase, out.receiver_z_m, 3.0 + zeros(size(cfg.receiver_x_m)));
verifySize(testCase, out.receiver_field, size(cfg.receiver_x_m));
end

function testInvalidReceiverCoordinatesFailBeforeAssembly(testCase)
cfg = local_config();
invalid = {struct('receiver_x_m',NaN,'receiver_z_m',3), ...
    struct('receiver_x_m',1i,'receiver_z_m',3), ...
    struct('receiver_x_m',0,'receiver_z_m',Inf), ...
    struct('receiver_x_m',[],'receiver_z_m',3)};
for ii = 1:numel(invalid)
    bad = cfg;
    bad.receiver_x_m = invalid{ii}.receiver_x_m;
    bad.receiver_z_m = invalid{ii}.receiver_z_m;
    did_fail = false;
    try
        solve_helmholtz_bie_halfplane_vertical(bad);
    catch exception
        did_fail = contains(exception.message,'Receiver coordinates must');
    end
    verifyTrue(testCase,did_fail);
end
end

function testOutOfWaterReceiverFailsBeforeAssembly(testCase)
cfg = local_config();
cfg.receiver_x_m = [0;0.1];
cfg.receiver_z_m = [3;0];
did_fail = false;
try
    solve_helmholtz_bie_halfplane_vertical(cfg);
catch exception
    did_fail = contains(exception.message,'Each receiver must lie in the water domain');
end
verifyTrue(testCase,did_fail);
end

function testMismatchedDepthVectorFailsBeforeAssembly(testCase)
cfg = local_config();
cfg.receiver_x_m = [-0.4; 0; 0.4];
cfg.receiver_z_m = [2.9; 3.1];

did_fail_for_depth_length = false;
try
    solve_helmholtz_bie_halfplane_vertical(cfg);
catch exception
    did_fail_for_depth_length = contains(exception.message, ...
        'receiver_z_m must be scalar or match receiver_x_m.');
end
verifyTrue(testCase, did_fail_for_depth_length);
end

function cfg = local_config()
frequency_hz = 1000;
k = 2*pi*frequency_hz/1500;
cfg = struct('frequency_hz', frequency_hz, 'c0_mps', 1500, ...
    'half_width_m', 5, 'points_per_wavelength', 4, ...
    'panel_order', 4, 'self_quadrature_order', 8, ...
    'halfplane_h_m', -2, 'window_plateau_ratio', 0.85, ...
    'eta_fn', @(x) zeros(size(x)), ...
    'eta_prime_fn', @(x) zeros(size(x)), ...
    'incident_fn', @(x, z) exp(1i*k.*z), ...
    'receiver_x_m', 0, 'receiver_z_m', 3);
end
