function receipt = audit_two_dimensional_closeout_s2c_b_phi_visual_semantics(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_PHI_VISUAL_SEMANTICS Compare saved FIG content.
% This supplements, but does not replace, the strict PNG pixel audit.
arguments
    run_dir (1,:) char
end
assert(isfolder(run_dir) && contains(run_dir,'20261023_s2c_b_phi_plus90_01'));
out=fullfile(run_dir,'audit','visual_semantics_phi_v1.json');
assert(~isfile(out),'Existing visual semantic audit receipt');
stems={ 'M05a_geometry_only','M05b_BIE_physical_response', ...
    'M05c_model_error_pair','M05d_line_center_nine','M05e_mirror_symmetry'};
result=cell(numel(stems),1);
for j=1:numel(stems)
    redraw=fullfile(run_dir,'audit','independent_redraw_phi_v1');
    if j==1
        redraw=fullfile(run_dir,'audit','independent_geometry_redraw_phi_v1');
    end
    original_fig=fullfile(run_dir,'figures_final_complete',[stems{j} '.fig']);
    redraw_fig=fullfile(redraw,[stems{j} '.fig']);
    assert(isfile(original_fig) && isfile(redraw_fig));
    f=openfig(original_fig,'invisible');
    g=openfig(redraw_fig,'invisible');
    cf=onCleanup(@()close(f)); cg=onCleanup(@()close(g));
    af=findall(f,'Type','axes'); ag=findall(g,'Type','axes');
    assert(numel(af)==numel(ag),'Axes count differs');
    nline=0; nimage=0;
    for k=1:numel(af)
        a=af(k); b=ag(k);
        assert(isequaln(a.XLim,b.XLim) && isequaln(a.YLim,b.YLim));
        assert(isequaln(a.Position,b.Position));
        assert(isequaln(a.Title.String,b.Title.String));
        assert(isequaln(a.XLabel.String,b.XLabel.String));
        assert(isequaln(a.YLabel.String,b.YLabel.String));
        x=findall(a,'Type','line'); y=findall(b,'Type','line');
        assert(numel(x)==numel(y),'Line count differs');
        nline=nline+numel(x);
        for p=1:numel(x)
            assert(isequaln(x(p).XData,y(p).XData));
            assert(isequaln(x(p).YData,y(p).YData));
            assert(isequaln(x(p).Color,y(p).Color));
            assert(isequaln(x(p).LineStyle,y(p).LineStyle));
            assert(isequaln(x(p).LineWidth,y(p).LineWidth));
        end
        x=findall(a,'Type','image'); y=findall(b,'Type','image');
        assert(numel(x)==numel(y),'Image count differs');
        nimage=nimage+numel(x);
        for p=1:numel(x)
            assert(isequaln(x(p).CData,y(p).CData));
            assert(isequaln(x(p).XData,y(p).XData));
            assert(isequaln(x(p).YData,y(p).YData));
        end
    end
    base_png=imread(fullfile(run_dir,'figures_final_complete',[stems{j} '.png']));
    new_png=imread(fullfile(redraw,[stems{j} '.png']));
    assert(isequal(size(base_png),size(new_png)));
    different_pixels=nnz(any(base_png~=new_png,3));
    result{j}=struct('stem',stems{j},'axes',numel(af),'lines',nline, ...
        'images',nimage,'fig_content_equal',true, ...
        'png_different_pixels',different_pixels,'png_total_pixels',size(base_png,1)*size(base_png,2));
    clear cf cg
end
receipt=struct('status','ALL_FIG_CONTENT_EQUAL_PNG_PIXEL_GATE_SEPARATE', ...
    'figures',{result},'strict_final_audit_passed',false, ...
    'note','Saved FIG axes, plotted data, nine-point image cells, labels and styles match. The original PNG pixel equality gate remains failed.');
fid=fopen(out,'w','n','UTF-8'); assert(fid>=0);
cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
end
