function fig = plot_bellhop_rays(ray_file, source_depth, receiver_depth, ...
    receiver_range_m, water_depth)
% plot_bellhop_rays - Plot ray paths from BELLHOP .ray file
%
% Returns empty if file not found

    fig = [];

    if ~isfile(ray_file)
        fprintf('Ray file not found: %s\n', ray_file);
        return;
    end

    fig = figure('Position', [250 250 1000 600]);
    hold on;

    % Read ray file
    fid = fopen(ray_file, 'r');

    % Skip header lines
    fgetl(fid);  % Title
    fgetl(fid);  % Frequency
    nsx_line = fgetl(fid);
    nsx = sscanf(nsx_line, '%d');
    fgetl(fid);  % nsy
    fgetl(fid);  % nsd

    ray_count = 0;
    max_range = 0;

    % Read rays
    for isx = 1:nsx
        nbeams_line = fgetl(fid);
        if ~ischar(nbeams_line), break; end
        nbeams = sscanf(nbeams_line, '%d');

        for ib = 1:nbeams
            npts_line = fgetl(fid);
            if ~ischar(npts_line), break; end
            npts = sscanf(npts_line, '%d');

            if npts > 0
                r = zeros(npts, 1);
                z = zeros(npts, 1);

                for ip = 1:npts
                    pt_line = fgetl(fid);
                    if ~ischar(pt_line), break; end
                    data = sscanf(pt_line, '%f');
                    if length(data) >= 2
                        r(ip) = data(1) * 1000;  % km to m
                        z(ip) = data(2);
                    end
                end

                if ray_count < 200
                    plot(r, z, 'b-', 'LineWidth', 0.5);
                    ray_count = ray_count + 1;
                    max_range = max(max_range, max(r));
                end
            end
        end
    end

    fclose(fid);

    % Plot boundaries
    rmax = max(receiver_range_m * 1.3, max_range * 1.1);
    plot([0 rmax], [0 0], 'k-', 'LineWidth', 2);  % Surface
    plot([0 rmax], [water_depth water_depth], 'k-', 'LineWidth', 2);  % Bottom

    % Plot source and receiver
    plot(0, source_depth, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Source');
    plot(receiver_range_m, receiver_depth, 'gs', 'MarkerSize', 12, ...
        'MarkerFaceColor', 'g', 'DisplayName', 'Receiver');

    set(gca, 'YDir', 'reverse');
    grid on;
    xlabel('Range (m)', 'FontSize', 12);
    ylabel('Depth (m)', 'FontSize', 12);
    title(sprintf('Ray Paths (%d rays)', ray_count), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best');
    xlim([0 rmax]);
    ylim([-10 water_depth + 10]);

    hold off;
end
