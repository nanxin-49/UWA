function H=deepsea_delay_separated_reconstruct(f_hz,path_response,delay_s,query_hz,meta)
%DEEPSEA_DELAY_SEPARATED_RECONSTRUCT Known paths/geometry, no path fitting.
arguments
    f_hz (:,1) double {mustBeFinite}
    path_response (:,:) double {mustBeFinite}
    delay_s (1,:) double {mustBeFinite}
    query_hz (:,1) double {mustBeFinite}
    meta (1,1) struct
end
assert(string(meta.physics_tier)=="analytic_control" && meta.known_path_components ...
    && string(meta.G1_target_PM)=="OPEN" && ~meta.WP2_formal, ...
    'Reconstruction:Scope','Requires analytic path components and frozen gates.');
assert(size(path_response,1)==numel(f_hz) && size(path_response,2)==numel(delay_s));
assert(all(diff(f_hz)>0) && min(query_hz)>=f_hz(1) && max(query_hz)<=f_hz(end));
slow=path_response.*exp(1i*2*pi*f_hz*delay_s);
coeff=interp1(f_hz,slow,query_hz,'linear');
H=sum(coeff.*exp(-1i*2*pi*query_hz*delay_s),2);
end
