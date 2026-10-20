function results = run_black_channel_basics(mode, saved_run_dir)
%RUN_BLACK_CHANNEL_BASICS BLACK-only download, official replay and basic analysis.
% setup_vertical_project; results = run_black_channel_basics();
% run_black_channel_basics('prepare') downloads/checks assets, never loads MAT.
% run_black_channel_basics('render', results.run_dir) redraws saved analysis.
if nargin < 1, mode = 'full'; end
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if strcmp(mode,'render')
    assert(nargin == 2,'black:RunDirectory','Saved run directory required.');
    loaded = load(fullfile(saved_run_dir,'analysis.mat'),'analysis');
    render_dir = fullfile(saved_run_dir,['render_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
    mkdir(render_dir);
    results = generate_black_channel_basics_visuals(loaded.analysis,render_dir);
    saved_status = jsondecode(fileread(fullfile(saved_run_dir,'status.json')));
    writetable(loaded.analysis.metrics,fullfile(render_dir,'black_metrics.csv'));
    write_black_channel_basics_report(loaded.analysis,saved_status,fullfile(render_dir,'black_report.md'));
    return
end
assert(any(strcmp(mode,{'full','prepare'})),'black:Mode','Use full, prepare or render.');
case_dir = fullfile(project_root,'results','validation','black_channel_basics');
run_dir = fullfile(case_dir,'runs',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
assert(~isfolder(run_dir),'black:ExistingRun','Refusing to overwrite a run.');
mkdir(run_dir);
results = struct('status','STARTED','run_dir',run_dir,'scope','BLACK only');
started = tic;
old_rng = rng;
rng_cleanup = onCleanup(@() rng(old_rng));
try
    results.matlab_version = version;
    results.current_directory = pwd;
    results.signal_toolbox = ver('signal');
    results.signal_license = license('test','Signal_Toolbox');
    assets = prepare_black_channel_assets(project_root);
    results.assets = assets.manifest;
    if strcmp(mode,'prepare')
        results.status = 'ASSETS_VERIFIED';
        local_json(fullfile(run_dir,'status.json'),results);
        return
    end
    assert(~isempty(results.signal_toolbox) && results.signal_license, ...
        'black:SignalDependency','Signal Processing Toolbox is required.');
    results.channel_inventory = black_inventory_mat(fullfile(assets.raw_dir,'black.mat'));
    results.noise_inventory = black_inventory_mat(fullfile(assets.raw_dir,'black_noise.mat'));
    results.status = 'INVENTORY_PASSED';
    local_json(fullfile(run_dir,'status.json'),results);
    channel = load(fullfile(assets.raw_dir,'black.mat'));
    noise = load(fullfile(assets.raw_dir,'black_noise.mat'));
    cfg = struct('seed',1994,'receivers',1:8,'representatives',[1 4 8], ...
        'fs_passband',100000,'symbol_rate',6250,'n_symbols',8192, ...
        'guard_s',0.02,'snr_db',20,'threshold_db',-25,'welch_window_s',4, ...
        'welch_overlap',0.5,'phase_block_s',0.1,'buffer_fraction',0.1);
    local_validate(channel,noise);
    results.synthetic = black_synthetic_control();
    original_path = path;
    path_cleanup = onCleanup(@() path(original_path));
    addpath(assets.vendor_dir,'-begin');
    assert(strcmp(which('replay'),fullfile(assets.vendor_dir,'replay.m')), ...
        'black:VendorPath','Pinned replay is shadowed.');
    fprintf('BLACK analysis: %d taps, %d receivers, %d snapshots\n',size(channel.h_hat));
    [metrics,pdp] = black_delay_metrics(channel.h_hat,channel.params.fs_delay);
    phase = black_phase_doppler(channel);
    M = size(channel.h_hat,2);
    if phase.available
        metrics.common_fd_median_hz = median(phase.frequency_hz,1).';
        metrics.common_fd_min_hz = min(phase.frequency_hz,[],1).';
        metrics.common_fd_max_hz = max(phase.frequency_hz,[],1).';
    else
        metrics.common_fd_median_hz = nan(M,1);
        metrics.common_fd_min_hz = nan(M,1);
        metrics.common_fd_max_hz = nan(M,1);
    end
    native_fs = channel.params.fs_time;
    K = size(channel.h_hat,1);
    buffer_left = cfg.buffer_fraction;
    buffer_right = cfg.buffer_fraction;
    if isfield(channel,'phi_hat') || isfield(channel,'f_resamp')
        n_phase = ceil(size(channel.h_hat,3)*channel.params.fs_delay/native_fs);
        drift_bounds = [0 0];
        for mm = cfg.representatives
            drift_phase = zeros(1,n_phase);
            if isfield(channel,'phi_hat'), drift_phase = channel.phi_hat(mm,:); end
            if isfield(channel,'f_resamp')
                drift_phase = drift_phase+(1/channel.f_resamp-1)*2*pi*channel.params.fc* ...
                    (0:n_phase-1)/channel.params.fs_delay;
            end
            drift_samples = drift_phase/(2*pi*channel.params.fc)*channel.params.fs_delay;
            drift_bounds = [min(drift_bounds(1),min(drift_samples)), ...
                max(drift_bounds(2),max(drift_samples))];
        end
        buffer_left = max(buffer_left,(ceil(drift_bounds(2))+2)/K);
        buffer_right = max(buffer_right,(ceil(-drift_bounds(1))+2)/K);
    end
    cfg.buffer_left = buffer_left;
    cfg.buffer_right = buffer_right;
    pad = ceil(K*buffer_left);
    % The official unpack inserts left padding; remove it from the plotted origin.
    h_reconstructed = unpack(native_fs,cfg.representatives,channel, ...
        buffer_left,buffer_right);
    tau_reconstructed = ((0:size(h_reconstructed,1)-1)-pad).'/channel.params.fs_delay;
    win_length = round(cfg.welch_window_s*native_fs);
    assert(size(channel.h_hat,3) >= win_length,'black:WelchDuration','Record shorter than four seconds.');
    residual_psd = [];
    for mm = 1:M
        tap = squeeze(channel.h_hat(metrics.main_tap_index(mm),mm,:));
        [power,f_doppler] = pwelch(tap,hann(win_length,'periodic'), ...
            floor(win_length/2),win_length,native_fs,'centered');
        residual_psd(:,mm) = power; %#ok<AGROW>
    end
    analysis = struct('cfg',cfg,'params',channel.params,'metrics',metrics,'pdp',pdp, ...
        'delay_s',(0:K-1).'/channel.params.fs_delay, ...
        'phase',phase,'residual_psd',residual_psd,'residual_frequency_hz',f_doppler, ...
        'h_reconstructed',h_reconstructed,'reconstructed_delay_s',tau_reconstructed, ...
        'reconstructed_time_s',(0:size(h_reconstructed,3)-1).'/native_fs, ...
        'duration_s',size(channel.h_hat,3)/native_fs,'manifest',assets.manifest);
    assert(all(isfinite(h_reconstructed(:))) && all(isfinite(residual_psd(:))), ...
        'black:AnalysisFinite','Nonfinite analysis array.');
    save(fullfile(run_dir,'analysis.mat'),'analysis','-v7.3');
    writetable(metrics,fullfile(run_dir,'black_metrics.csv'));
    results.status = 'ANALYSIS_SAVED';
    local_json(fullfile(run_dir,'status.json'),results);
    fprintf('Official BLACK replay, all eight receivers ...\n');
    rng(cfg.seed,'twister');
    symbols = 2*randi([0 1],cfg.n_symbols,1)-1;
    baseband = resample(symbols,cfg.fs_passband/cfg.symbol_rate,1);
    time_input = (0:numel(baseband)-1).'/cfg.fs_passband;
    passband = real(baseband.*exp(2j*pi*channel.params.fc*time_input));
    guard = zeros(round(cfg.guard_s*cfg.fs_passband),1);
    input = [guard;passband;guard];
    start = 1+round(0.25*analysis.duration_s*channel.params.fs_delay);
    y = replay(input,cfg.fs_passband,cfg.receivers,channel,start);
    y_repeat = replay(input,cfg.fs_passband,cfg.receivers,channel,start);
    repeat_relative_error = norm(y_repeat(:)-y(:))/max(norm(y(:)),realmin);
    assert(repeat_relative_error <= 1e-12,'black:Repeat','Repeat replay mismatch.');
    clear y_repeat
    rng(cfg.seed,'twister');
    w = noisegen(size(y),cfg.fs_passband,cfg.receivers,noise);
    sigma = sqrt(mean(y(:).^2)/(10^(cfg.snr_db/10)*mean(w(:).^2)));
    noise_added = sigma*w;
    received = y+noise_added;
    measured_snr_db = 10*log10(mean(y(:).^2)/mean(noise_added(:).^2));
    assert(abs(measured_snr_db-cfg.snr_db)<1e-10,'black:SNR','SNR scaling failed.');
    assert(all(isfinite(y(:))) && all(isfinite(received(:))), 'black:ReplayFinite','Nonfinite replay.');
    arrival = black_check_replay_arrival(y,baseband,channel,cfg,start);
    save(fullfile(run_dir,'replay.mat'),'cfg','input','symbols','baseband','y','received', ...
        'noise_added','sigma','start','repeat_relative_error','measured_snr_db','arrival','-v7.3');
    results.repeat_relative_error = repeat_relative_error;
    results.measured_snr_db = measured_snr_db;
    results.arrival = arrival;
    results.cfg = cfg;
    results.figures = generate_black_channel_basics_visuals(analysis,run_dir);
    if arrival.passed, results.status = 'PASS';
    else, results.status = 'REPLAY_ARRIVAL_FAILED'; end
    results.elapsed_s = toc(started);
    results.directory_unchanged = strcmp(results.current_directory,pwd);
    source_names = {mfilename('fullpath'),'black_delay_metrics','black_phase_doppler', ...
        'black_check_replay_arrival','black_inventory_mat','black_synthetic_control', ...
        'prepare_black_channel_assets','black_file_hash', ...
        'generate_black_channel_basics_visuals','write_black_channel_basics_report', ...
        fullfile(project_root,'scripts','reporting','templates','black_channel_basics_report_zh.md')};
    workflow_dir = fullfile(run_dir,'workflow_source');
    mkdir(workflow_dir);
    results.workflow_code = repmat(struct('name','','sha256',''),1,numel(source_names));
    for ii=1:numel(source_names)
        if ii==1
            source=[source_names{ii} '.m'];
        elseif isfile(source_names{ii})
            source=source_names{ii};
        else
            source=which(source_names{ii});
        end
        [~,name,ext] = fileparts(source);
        copyfile(source,fullfile(workflow_dir,[name ext]));
        results.workflow_code(ii) = struct('name',[name ext],'sha256',black_file_hash(source));
    end
    local_json(fullfile(run_dir,'status.json'),results);
    write_black_channel_basics_report(analysis,results,fullfile(run_dir,'black_report.md'));
    write_black_channel_basics_report(analysis,results, ...
        fullfile(project_root,'reports','black_channel_basics_report.md'));
    fprintf('BLACK result: %s; run directory: %s\n',results.status,run_dir);
catch ME
    results.status = 'FAILED';
    results.error_identifier = ME.identifier;
    results.error_message = ME.message;
    results.error_stack = ME.stack;
    results.elapsed_s = toc(started);
    local_json(fullfile(run_dir,'status.json'),results);
    rethrow(ME);
end
end

function local_validate(channel,noise)
assert(isfield(channel,'params') && isfield(channel,'h_hat') && channel.version>=1, ...
    'black:Format','Required channel fields missing.');
assert(size(channel.h_hat,2)==8 && abs(channel.params.fc-18000)<1, ...
    'black:Identity','Not the expected BLACK channel.');
assert(all(isfinite(channel.h_hat(:))) && all(isfinite(noise.beta(:))), ...
    'black:FiniteInput','Nonfinite raw data.');
assert(~(isfield(channel,'phi_hat') && isfield(channel,'theta_hat')), ...
    'black:TrackingAmbiguity','Official unpack/replay differ when both phase fields coexist.');
for field = {'phi_hat','theta_hat'}
    if isfield(channel,field{1})
        phase = channel.(field{1});
        assert(size(phase,1)==8 && all(isfinite(phase(:))), 'black:Phase','Invalid phase array.');
        discrepancy = abs(size(phase,2)/channel.params.fs_delay- ...
            size(channel.h_hat,3)/channel.params.fs_time);
        assert(discrepancy <= 1/channel.params.fs_time,'black:Duration','Inconsistent durations.');
    end
end
assert(size(noise.beta,1)==8 && size(noise.beta,2)==8 && noise.alpha==2, ...
    'black:Noise','Unexpected BLACK noise model.');
assert(abs(noise.fc-channel.params.fc)<1,'black:NoiseCarrier','Noise/channel carriers differ.');
end

function local_json(filename,value)
fid = fopen(filename,'w','n','UTF-8');
assert(fid>=0,'black:WriteStatus','Cannot write status.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
