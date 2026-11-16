function [meta,kind]=write_bellhop_rotated_REAL8_C_inputs(root,cfg)
%WRITE_BELLHOP_ROTATED_REAL8_B_INPUTS Unchanged validated writers, no label adjustment.
if isfield(cfg,'reflection_order')
    kind='reflection'; meta=write_bellhop_rotated_REAL8_A_inputs(root,cfg);
elseif strcmp(cfg.chart_kind,'horizontal')
    kind='horizontal'; meta=write_bellhop_rotated_horizontal_q(root,cfg,cfg.horizontal_side); meta.q_receiver_m=cfg.receiver_q_m;
else
    kind='direct'; meta=write_bellhop_rotated_munk_q(root,cfg);
    fid=fopen([root '.iwcov'],'w'); assert(fid>=0); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
    write_bellhop_rotated_field_options(root,cfg.receiver_range_m,meta.q_receiver_m,point3d=true,fast_q=true,native_range_metadata_repair=false);
end
end