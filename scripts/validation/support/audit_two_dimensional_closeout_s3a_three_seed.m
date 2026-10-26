function result = audit_two_dimensional_closeout_s3a_three_seed(root,dest)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S3A_THREE_SEED Independent saved summary audit.
arguments
    root (1,:) char
    dest (1,:) char
end
parent=fullfile(root,'results','comparisons','two_dimensional_channel_closeout');
runs={'20261025_s3a_seed260001_01','20261025_s3a_seed260002_01', ...
    '20261026_s3a_seed260003_01'};
seeds=[260001 260002 260003];
table_path=fullfile(dest,'three_seed_comparison.csv');
geo_path=fullfile(dest,'three_seed_geometry.csv');
summary=readtable(table_path,TextType='string');
geometry=readtable(geo_path,TextType='string');
assert(height(summary)==6 && height(geometry)==3);
max_difference=0;
for k=1:3
    run=fullfile(parent,runs{k});
    input=jsondecode(fileread(fullfile(run,'audit','input_audit.json')));
    assert(input.seeds(3).seed==260003 && input.seeds(3).same_parent_max_abs==0 ...
        && input.seeds(3).source_coeff_max_abs==0);
    m=readtable(fullfile(run,'analysis_v2','model_metrics.csv'),TextType='string');
    p=readtable(fullfile(run,'analysis_v2','BIE_physical_response.csv'),TextType='string');
    for j=1:2
        sid=["nine_9" "line_411"]; sid=sid(j);
        row=summary(summary.seed==seeds(k) & summary.support_id==sid,:);
        pe=m(m.method=="PE" & m.support_id==sid,:);
        fk=m(m.method=="FullFK" & m.support_id==sid,:);
        pb=p(p.support_id==sid,:);
        assert(height(row)==1 && height(pe)==1 && height(fk)==1 && height(pb)==1);
        ref=[pb.BIE_physical_response_L2,pb.BIE_level_resolution_relative, ...
            pe.complex_L2,pe.amplitude_L2,pe.phase_L2,pe.phase_RMS_rad, ...
            pe.total_empirical_resolution,fk.complex_L2,fk.amplitude_L2, ...
            fk.phase_L2,fk.phase_RMS_rad,fk.total_empirical_resolution];
        current=[row.BIE_physical_response_L2,row.BIE_reference_resolution, ...
            row.PE_complex_L2,row.PE_amplitude_L2,row.PE_phase_L2, ...
            row.PE_phase_RMS_rad,row.PE_numeric_resolution, ...
            row.FullFK_complex_L2,row.FullFK_amplitude_L2,row.FullFK_phase_L2, ...
            row.FullFK_phase_RMS_rad,row.FullFK_numeric_resolution];
        max_difference=max(max_difference,max(abs(ref-current)));
        assert(row.n==numel(find(strcmp(sid,["nine_9" "line_411"])))+ ...
            (sid=="nine_9")*8+(sid=="line_411")*410);
        assert(abs(row.PE_complex_L2^2-row.PE_amplitude_L2^2- ...
            row.PE_phase_L2^2)<1e-12);
        assert(abs(row.FullFK_complex_L2^2-row.FullFK_amplitude_L2^2- ...
            row.FullFK_phase_L2^2)<1e-12);
        assert(row.PE_complex_L2>5*row.PE_numeric_resolution && ...
            row.FullFK_complex_L2>5*row.FullFK_numeric_resolution);
    end
    assert(geometry.seed(k)==seeds(k) && geometry.embedded_min_height_m(k)>-2);
end
assert(max_difference<1e-12);
fig=fullfile(dest,'F08_S3A_three_seed_summary.fig');
png=fullfile(dest,'F08_S3A_three_seed_summary.png');
h=openfig(fig,'invisible');
axes=findall(h,'Type','axes');
assert(numel(axes)==3 && numel(findall(h,'Type','line'))==6);
redraw=fullfile(dest,'F08_S3A_three_seed_summary_fig_reexport.png');
exportgraphics(h,redraw,'Resolution',160); close(h);
assert(isequal(imread(png),imread(redraw)));
result=struct('status','S3A_THREE_SEED_INDEPENDENT_AUDIT_PASS', ...
    'seed_count',3,'support_rows',6,'max_saved_CSV_difference',max_difference, ...
    'F08_FIG_PNG_pixel_difference',0, ...
    'qualification','three fixed PM realizations only');
file=fullfile(dest,'independent_audit.json');
fid=fopen(file,'w'); assert(fid>0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(result,PrettyPrint=true));
end
