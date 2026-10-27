function tests = test_two_dimensional_s5_bie_receiver_coordinates
%TEST_TWO_DIMENSIONAL_S5_BIE_RECEIVER_COORDINATES Receiver evaluation contract.
tests = functiontests(localfunctions);
end

function testPairedCoordinatesMatchSeparateEvaluations(testCase)
cfg = local_config();
cfg.receiver_x_m = [-0.4; 0; 0.4];
cfg.receiver_z_m = [2.9; 3.0; 3.1];
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
end

function testScalarDepthStillExpandsAcrossReceivers(testCase)
cfg = local_config();
cfg.receiver_x_m = [-0.4; 0; 0.4];
cfg.receiver_z_m = 3.0;
out = solve_helmholtz_bie_halfplane_vertical(cfg);

verifyEqual(testCase, out.receiver_z_m, 3.0 + zeros(size(cfg.receiver_x_m)));
verifySize(testCase, out.receiver_field, size(cfg.receiver_x_m));
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
