function result = run_bellhop_flat_surface_arrivals_vertical(cfg)
%RUN_BELLHOP_FLAT_SURFACE_ARRIVALS_VERTICAL Execute one flat-surface case.
% Validation-only helper. The caller owns the output directory and metadata.

arguments
    cfg (1,1) struct
end
required={'bellhop_exe','case_root','frequency_hz','water_depth_m','c0_mps', ...
    'z_tx_m','z_rx_m','offset_m','beam_count','angle_min_deg','angle_max_deg', ...
    'arrival_cluster_tolerance_ms'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii})
        error('run_bellhop_flat_surface_arrivals_vertical:MissingField', ...
            'cfg.%s is required.',required{ii});
    end
end
validateattributes(cfg.offset_m,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.beam_count,{'numeric'},{'scalar','integer','>=',2});
[~,~,case_extension]=fileparts(cfg.case_root);
if ~isempty(case_extension)
    error('run_bellhop_flat_surface_arrivals_vertical:InvalidCaseRoot', ...
        'cfg.case_root must not contain a filename extension: %s',cfg.case_root);
end
if exist(cfg.bellhop_exe,'file')~=2
    error('run_bellhop_flat_surface_arrivals_vertical:MissingExecutable', ...
        'Bellhop executable does not exist: %s',cfg.bellhop_exe);
end

out_dir=fileparts(cfg.case_root);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
env_file=[cfg.case_root '.env'];
arr_file=[cfg.case_root '.arr'];
prt_file=[cfg.case_root '.prt'];
local_write_env(env_file,cfg);
if exist(arr_file,'file')==2, delete(arr_file); end

old=pwd; cleanup=onCleanup(@()cd(old)); cd(out_dir);
[~,name]=fileparts(cfg.case_root);
[status,command_output]=system(sprintf('"%s" "%s"',cfg.bellhop_exe,name));
clear cleanup
if status~=0 || exist(arr_file,'file')~=2
    error('run_bellhop_flat_surface_arrivals_vertical:ExecutionFailed', ...
        'Bellhop failed for %s: %s',name,command_output);
end

raw=local_read_arrivals(arr_file);
geometry=struct( ...
    'direct_length_m',hypot(cfg.offset_m,cfg.z_tx_m-cfg.z_rx_m), ...
    'surface_length_m',hypot(cfg.offset_m,cfg.z_tx_m+cfg.z_rx_m));
geometry.direct_time_s=geometry.direct_length_m/cfg.c0_mps;
geometry.surface_time_s=geometry.surface_length_m/cfg.c0_mps;
paths=local_classify(raw,geometry,cfg.arrival_cluster_tolerance_ms);
result=struct('config',cfg,'geometry',geometry,'raw_arrivals',raw, ...
    'paths',paths,'command_output',command_output,'files',struct( ...
    'env',env_file,'arr',arr_file,'prt',prt_file));
end

function local_write_env(file,cfg)
fid=fopen(file,'w');
if fid<0, error('Cannot create Bellhop environment: %s',file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'''PE-Bellhop current flat-surface validation''\n%.12g\n1\n''CVW''\n', ...
    cfg.frequency_hz);
fprintf(fid,'2 0 %.12g\n0 %.12g /\n%.12g %.12g /\n', ...
    cfg.water_depth_m,cfg.c0_mps,cfg.water_depth_m,cfg.c0_mps);
fprintf(fid,'''A'' 0\n%.12g 1800 0 2 0 /\n',cfg.water_depth_m);
fprintf(fid,'1\n%.12g /\n1\n%.12g /\n1\n%.12g /\n', ...
    cfg.z_tx_m,cfg.z_rx_m,cfg.offset_m/1000);
fprintf(fid,'''A''\n%d\n%.12g %.12g /\n0 %.12g %.12g\n', ...
    cfg.beam_count,cfg.angle_min_deg,cfg.angle_max_deg, ...
    1.1*cfg.water_depth_m,max(2*cfg.offset_m/1000,1e-6));
clear cleanup
end

function t=local_read_arrivals(file)
fid=fopen(file,'r');
if fid<0, error('Cannot open Bellhop arrivals: %s',file); end
cleanup=onCleanup(@()fclose(fid));
f=fscanf(fid,'%f',1); ns=fscanf(fid,'%d',1); nd=fscanf(fid,'%d',1); nr=fscanf(fid,'%d',1);
sd=fscanf(fid,'%f',ns); rd=fscanf(fid,'%f',nd); rr=1000*fscanf(fid,'%f',nr);
if ns~=1 || nd~=1 || nr~=1
    error('Expected one source depth, receiver depth, and range.');
end
fscanf(fid,'%d',1); n=fscanf(fid,'%d',1); values=fscanf(fid,'%f',[8,n]);
if size(values,2)~=n, error('Incomplete Bellhop arrivals file: %s',file); end
amplitude_complex=values(1,:).'.*exp(1i*deg2rad(values(2,:).'));
t=table((1:n).',repmat(f,n,1),repmat(sd,n,1),repmat(rd,n,1), ...
    repmat(rr,n,1),amplitude_complex,abs(amplitude_complex),values(2,:).', ...
    values(3,:).',values(4,:).',values(5,:).',values(6,:).', ...
    values(7,:).',values(8,:).','VariableNames', ...
    {'arrival_index','frequency_hz','source_depth_m','receiver_depth_m', ...
    'receiver_range_m','amplitude_complex','amplitude','phase_deg','delay_s', ...
    'delay_imag_s','source_angle_deg','receiver_angle_deg', ...
    'top_bounce_count','bottom_bounce_count'});
clear cleanup
end

function paths=local_classify(raw,geometry,tolerance_ms)
names=["direct";"surface_reflection"];
top_counts=[0;1]; targets=[geometry.direct_time_s;geometry.surface_time_s];
template=table(names,nan(2,1),complex(nan(2,1)),nan(2,1),nan(2,1), ...
    nan(2,1),zeros(2,1),top_counts,zeros(2,1),'VariableNames', ...
    {'path_name','delay_s','amplitude_complex','amplitude','source_angle_deg', ...
    'receiver_angle_deg','cluster_size','top_bounce_count','bottom_bounce_count'});
paths=template;
for pp=1:2
    q=sortrows(raw(raw.bottom_bounce_count==0 & ...
        raw.top_bounce_count==top_counts(pp),:),'delay_s');
    if isempty(q)
        error('Bellhop did not return the %s path.',names(pp));
    end
    breaks=[true;diff(q.delay_s)>tolerance_ms/1000];
    ids=cumsum(breaks); unique_ids=unique(ids);
    centers=arrayfun(@(id)mean(q.delay_s(ids==id)),unique_ids);
    [~,best]=min(abs(centers-targets(pp))); mask=ids==unique_ids(best);
    weights=max(q.amplitude(mask),realmin);
    paths.delay_s(pp)=sum(weights.*q.delay_s(mask))/sum(weights);
    paths.amplitude_complex(pp)=sum(q.amplitude_complex(mask));
    paths.amplitude(pp)=abs(paths.amplitude_complex(pp));
    paths.source_angle_deg(pp)=sum(weights.*q.source_angle_deg(mask))/sum(weights);
    paths.receiver_angle_deg(pp)=sum(weights.*q.receiver_angle_deg(mask))/sum(weights);
    paths.cluster_size(pp)=sum(mask);
end
end
