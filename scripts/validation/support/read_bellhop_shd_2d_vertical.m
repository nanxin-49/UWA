function shade = read_bellhop_shd_2d_vertical(file)
%READ_BELLHOP_SHD_2D_VERTICAL Read a rectilinear 2-D Bellhop shade file.
fid=fopen(file,'rb');
if fid<0, error('Cannot open Bellhop shade file: %s',file); end
cleanup=onCleanup(@()fclose(fid));
record_words=fread(fid,1,'int32');
if isempty(record_words) || record_words<=0, error('Invalid shade record length.'); end
record_bytes=4*record_words;
title_text=deblank(fread(fid,80,'*char').');
fseek(fid,record_bytes,'bof'); plot_type=deblank(fread(fid,10,'*char').');
fseek(fid,2*record_bytes,'bof');
frequency_hz=fread(fid,1,'float32');
counts=fread(fid,6,'int32').'; attenuation=fread(fid,1,'float32');
if numel(counts)~=6 || counts(1)~=1 || counts(4)~=1
    error('Expected single-bearing, single-source-depth 2-D shade data.');
end
n_theta=counts(1); n_sx=counts(2); n_sy=counts(3);
n_sd=counts(4); n_rd=counts(5); n_rr=counts(6);
fseek(fid,3*record_bytes,'bof'); theta_deg=fread(fid,n_theta,'float32');
fseek(fid,4*record_bytes,'bof'); source_x_m=fread(fid,n_sx,'float32');
fseek(fid,5*record_bytes,'bof'); source_y_m=fread(fid,n_sy,'float32');
fseek(fid,6*record_bytes,'bof'); source_depth_m=fread(fid,n_sd,'float32');
fseek(fid,7*record_bytes,'bof'); receiver_depth_m=fread(fid,n_rd,'float32');
fseek(fid,9*record_bytes,'bof'); receiver_range_m=fread(fid,n_rr,'float32');
pressure=complex(zeros(n_rd,n_rr));
for dd=1:n_rd
    fseek(fid,(9+dd)*record_bytes,'bof'); raw=fread(fid,2*n_rr,'float32');
    if numel(raw)~=2*n_rr, error('Incomplete shade pressure record %d.',dd); end
    pressure(dd,:)=raw(1:2:end)+1i*raw(2:2:end);
end
shade=struct('title',title_text,'plot_type',plot_type, ...
    'frequency_hz',frequency_hz,'attenuation',attenuation, ...
    'theta_deg',theta_deg,'source_x_m',source_x_m,'source_y_m',source_y_m, ...
    'source_depth_m',source_depth_m,'receiver_depth_m',receiver_depth_m, ...
    'receiver_range_m',receiver_range_m,'pressure',pressure, ...
    'counts',counts,'record_bytes',record_bytes,'file',file);
clear cleanup
end
