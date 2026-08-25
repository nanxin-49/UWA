function shade = read_bellhop_shd_unfolded_vertical(file)
%READ_BELLHOP_SHD_UNFOLDED_VERTICAL Read Bellhop SHD with optional .sbp.
% Bellhop versions used by this project shift the rectilinear range/data
% records by one fixed record when a source beam pattern is enabled.  This
% reader identifies the range record from its monotone positive values and
% then reads the depth blocks relative to it.
fid=fopen(file,'rb'); if fid<0, error('Cannot open Bellhop shade file: %s',file); end
cleanup=onCleanup(@()fclose(fid)); record_words=fread(fid,1,'int32');
if isempty(record_words) || record_words<=0, error('Invalid shade record length.'); end
rb=4*record_words; title_text=deblank(fread(fid,80,'*char').');
fseek(fid,2*rb,'bof'); frequency_hz=fread(fid,1,'float32'); counts=fread(fid,6,'int32').'; attenuation=fread(fid,1,'float32');
if numel(counts)~=6 || counts(1)~=1 || counts(4)~=1, error('Expected single-bearing, single-source-depth SHD.'); end
n_theta=counts(1); n_sx=counts(2); n_sy=counts(3); n_sd=counts(4); n_rd=counts(5); n_rr=counts(6);
fseek(fid,3*rb,'bof'); theta_deg=fread(fid,n_theta,'float32');
fseek(fid,4*rb,'bof'); source_x_m=fread(fid,n_sx,'float32');
fseek(fid,5*rb,'bof'); source_y_m=fread(fid,n_sy,'float32');
fseek(fid,6*rb,'bof'); source_depth_m=fread(fid,n_sd,'float32');
fseek(fid,7*rb,'bof'); receiver_depth_m=fread(fid,n_rd,'float32');
range_record=local_find_range_record(fid,rb,n_rr); fseek(fid,range_record*rb,'bof'); receiver_range_m=fread(fid,n_rr,'float32');
pressure=complex(zeros(n_rd,n_rr));
for dd=1:n_rd
    fseek(fid,(range_record+dd)*rb,'bof'); raw=fread(fid,2*n_rr,'float32');
    if numel(raw)~=2*n_rr, error('Incomplete unfolded SHD pressure record %d.',dd); end
    pressure(dd,:)=raw(1:2:end)+1i*raw(2:2:end);
end
shade=struct('title',title_text,'plot_type','rectilin','frequency_hz',frequency_hz, ...
    'attenuation',attenuation,'theta_deg',theta_deg,'source_x_m',source_x_m, ...
    'source_y_m',source_y_m,'source_depth_m',source_depth_m,'receiver_depth_m',receiver_depth_m, ...
    'receiver_range_m',receiver_range_m,'pressure',pressure,'counts',counts, ...
    'record_bytes',rb,'range_record',range_record,'file',file);
clear cleanup
end

function index=local_find_range_record(fid,rb,n_rr)
candidates=[8 9 10];
for ii=1:numel(candidates)
    fseek(fid,candidates(ii)*rb,'bof'); q=fread(fid,n_rr,'float32');
    if numel(q)==n_rr && all(isfinite(q)) && all(q>0) && all(diff(q)>=0) && q(end)>q(1)
        index=candidates(ii); return
    end
end
error('Unable to locate monotone positive receiver-range record.');
end
