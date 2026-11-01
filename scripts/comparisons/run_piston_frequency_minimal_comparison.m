function result = run_piston_frequency_minimal_comparison(run_id)
%RUN_PISTON_FREQUENCY_MINIMAL_COMPARISON Analytic piston frequency sweep.
if nargin < 1 || isempty(run_id)
    run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));
end
run_id = char(run_id);
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')));
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(root,'results','comparisons','source_model_minimal',run_id);
assert(~exist(out,'dir'),'Run directory already exists.');
mkdir(out);
cfg = struct('radius_m',0.05,'c_mps',1500,'frequencies_hz',[4000 8000 12000], ...
    'theta_marker_deg',7.6756);
theta_deg = (0:0.01:10).';
theta_rad = theta_deg*pi/180;
D = zeros(numel(theta_deg),numel(cfg.frequencies_hz));
marker_db = zeros(numel(cfg.frequencies_hz),1);
for ii = 1:numel(cfg.frequencies_hz)
    ka = 2*pi*cfg.frequencies_hz(ii)*cfg.radius_m/cfg.c_mps;
    u = ka*sin(theta_rad);
    D(:,ii) = local_piston(u);
    marker_db(ii) = 20*log10(local_piston(ka*sin(cfg.theta_marker_deg*pi/180)));
end
assert(all(isfinite(D(:))) && all(D(1,:)==1),'Piston axis check failed.');
assert(all(diff(marker_db)<0),'Expected monotonic narrowing at the marked angle.');
T = table(cfg.frequencies_hz(:)/1000,marker_db, ...
    'VariableNames',{'frequency_khz','level_at_7p6756_deg_db'});
curve_table = table(theta_deg,D(:,1),D(:,2),D(:,3), ...
    'VariableNames',{'theta_deg','D_4kHz','D_8kHz','D_12kHz'});
writetable(T,fullfile(out,'frequency_summary.csv'));
writetable(curve_table,fullfile(out,'direction_curves.csv'));
fig = figure('Visible','off','Color','w','Position',[100 100 920 550]);
ax = axes(fig); hold(ax,'on');
styles = {'-','--','-.'};
for ii = 1:3
    plot(ax,theta_deg,20*log10(D(:,ii)),styles{ii},'LineWidth',2, ...
        'DisplayName',sprintf('%g kHz',cfg.frequencies_hz(ii)/1000));
end
line_marker = xline(ax,cfg.theta_marker_deg,':','7.6756 deg','LineWidth',1.5);
line_marker.HandleVisibility = 'off';
xlabel(ax,'Angle from upward axis (deg)');
ylabel(ax,'Axis-normalized pressure level (dB)');
title(ax,'Ideal 5 cm circular piston: frequency dependence near the axis');
xlim(ax,[0 10]); ylim(ax,[-0.25 0.01]); grid(ax,'on');
legend(ax,'Location','southwest');
drawnow;
savefig(fig,fullfile(out,'piston_frequency_pattern.fig'));
exportgraphics(fig,fullfile(out,'piston_frequency_pattern.png'),'Resolution',180);
close(fig);
result = struct('run_id',run_id,'config',cfg,'summary',T,'curve',curve_table, ...
    'failed_checks',0,'output_dir',out);
save(fullfile(out,'piston_frequency_data.mat'),'result','-v7');
disp(T);
end

function D = local_piston(u)
D = ones(size(u));
nonzero = u~=0;
D(nonzero) = 2*besselj(1,u(nonzero))./u(nonzero);
end
