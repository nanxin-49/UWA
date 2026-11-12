function c=munk_profile_sound_speed(z,p)
%MUNK_PROFILE_SOUND_SPEED Existing Fig12b Munk formula, explicit parameters.
arguments
    z double {mustBeFinite}
    p (1,1) struct
end
a=2*(z-p.axis_depth_m)/p.scale_m;
c=p.c0_mps*(1+p.epsilon*(a-1+exp(-a)));
end
