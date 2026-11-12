function tau=munk_vertical_travel_time_reference(z,source_z,p)
%MUNK_VERTICAL_TRAVEL_TIME_REFERENCE Independent vertical ray time integral.
% Added after current-tree search failed and user explicitly permitted it.
% This is not the full multi-path field's group delay or a wave-solver PASS.
arguments
    z double {mustBeFinite,mustBeNonnegative}
    source_z (1,1) double {mustBeFinite,mustBePositive}
    p (1,1) struct
end
tau=zeros(size(z));
for j=1:numel(z)
    tau(j)=abs(integral(@(v)1./munk_profile_sound_speed(v,p),source_z,z(j), ...
        'RelTol',1e-12,'AbsTol',1e-13));
end
end
