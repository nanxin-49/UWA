function [arrivals, summary] = read_bellhop_arrivals(arr_file)
% read_bellhop_arrivals - Read BELLHOP .arr file
%
% Returns:
%   arrivals - struct with fields: amplitude, delay, src_angle, rcv_angle,
%              num_top_bounces, num_bot_bounces
%   summary - struct with statistics

    arrivals = struct();
    summary = struct();
    summary.num_arrivals = 0;

    if ~isfile(arr_file)
        fprintf('Arrivals file not found: %s\n', arr_file);
        return;
    end

    fid = fopen(arr_file, 'r');

    % Read header
    fgetl(fid);  % Title
    freq_line = fgetl(fid);
    freq = sscanf(freq_line, '%f');

    nsd_line = fgetl(fid);
    nsd = sscanf(nsd_line, '%d');

    nrd_line = fgetl(fid);
    nrd = sscanf(nrd_line, '%d');

    nrr_line = fgetl(fid);
    nrr = sscanf(nrr_line, '%d');

    % Initialize arrays
    amplitude = [];
    delay = [];
    src_angle = [];
    rcv_angle = [];
    num_top_bounces = [];
    num_bot_bounces = [];

    % Read arrivals
    for isd = 1:nsd
        for ird = 1:nrd
            for irr = 1:nrr
                narr_line = fgetl(fid);
                if ~ischar(narr_line), break; end
                narr = sscanf(narr_line, '%d');
                if isempty(narr), narr = 0; end

                for ia = 1:narr
                    arr_line = fgetl(fid);
                    if ~ischar(arr_line), break; end

                    data = sscanf(arr_line, '%f');
                    if length(data) >= 8
                        amplitude = [amplitude; data(1) + 1i*data(2)];
                        delay = [delay; data(4)];
                        src_angle = [src_angle; data(5)];
                        rcv_angle = [rcv_angle; data(6)];
                        num_top_bounces = [num_top_bounces; data(7)];
                        num_bot_bounces = [num_bot_bounces; data(8)];
                    end
                end
            end
        end
    end

    fclose(fid);

    % Store results
    arrivals.amplitude = amplitude;
    arrivals.delay = delay;
    arrivals.src_angle = src_angle;
    arrivals.rcv_angle = rcv_angle;
    arrivals.num_top_bounces = num_top_bounces;
    arrivals.num_bot_bounces = num_bot_bounces;

    % Compute summary
    num_arr = length(amplitude);
    summary.num_arrivals = num_arr;

    if num_arr > 0
        summary.first_delay = min(delay);
        summary.last_delay = max(delay);
        summary.delay_spread = summary.last_delay - summary.first_delay;

        [max_amp, max_idx] = max(abs(amplitude));
        summary.max_amplitude = max_amp;
        summary.max_delay = delay(max_idx);
    else
        summary.first_delay = NaN;
        summary.last_delay = NaN;
        summary.delay_spread = NaN;
        summary.max_amplitude = NaN;
        summary.max_delay = NaN;
    end
end
