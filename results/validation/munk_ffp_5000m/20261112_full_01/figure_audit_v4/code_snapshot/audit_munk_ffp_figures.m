function receipt=audit_munk_ffp_figures(run_id,version_suffix,audit_id)
%AUDIT_MUNK_FFP_FIGURES Data/coordinate audit plus FIG reopen and PNG roundtrip.
arguments
    run_id (1,:) char = '20261112_full_01'
    version_suffix (1,:) char = '_v3'
    audit_id (1,:) char = 'figure_audit_v3'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(isempty(regexp([run_id version_suffix audit_id],'[^A-Za-z0-9_-]','once')));
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id);
v=fullfile(root,'results','visualization','munk_ffp_5000m',[run_id version_suffix]);
a=fullfile(d,audit_id); assert(~isfolder(a)); mkdir(a); t=tic;
receipt=struct('stage','figure_audit','status','RUNNING','elapsed_s',0,'exit_status',NaN,'failed_checks',0);
try
    s=load(fullfile(d,'J07.mat'),'output'); o=s.output;
    stems={'F01_full_field_vs_paper','F02_axis_magnitude_phase','F03_near_bottom_resolution','F04_uniform_reference_errors','F05_convergence_qualification'};
    f=openfig(fullfile(v,[stems{1} '.fig']),'new','invisible'); images=findobj(f,'Type','image');
    globalz=(0:5:5000)'; [present,ind]=ismember(globalz,o.depth_m); assert(all(present));
    expected=-20*log10(abs(o.pressure_1m(ind,:))); expected(~isfinite(expected))=NaN; found=false;
    for j=1:numel(images)
        if ndims(images(j).CData)==2
            assert(isequal(size(images(j).CData),size(expected)), ...
                'munk_ffp:NonuniformImageMapping','F01 must use actual uniform global rows (1001x201), not combined nonuniform rows');
            assert(isequaln(images(j).CData,expected));
            assert(abs(images(j).YData(1)-0)<1e-12 && abs(images(j).YData(end)-5000)<1e-12);
            found=true;
        end
    end
    close(f); assert(found);
    f=openfig(fullfile(v,[stems{2} '.fig']),'new','invisible'); ls=findobj(f,'Type','line'); axis=o.range_m==0;
    amplitude=abs(o.pressure_1m(:,axis)); magnitude_found=false;
    for j=1:numel(ls)
        if numel(ls(j).XData)==numel(amplitude) && isequaln(ls(j).YData(:),o.depth_m(:)) && isequaln(ls(j).XData(:),amplitude(:))
            magnitude_found=true;
        end
    end
    close(f); assert(magnitude_found,'Axial actual data not found in F02');
    f=openfig(fullfile(v,[stems{3} '.fig']),'new','invisible'); images=findobj(f,'Type','image');
    actual=0;
    for j=1:numel(images)
        yy=images(j).YData; rows=o.depth_m>=min(yy) & o.depth_m<=max(yy);
        expected=-20*log10(abs(o.pressure_1m(rows,:))); expected(~isfinite(expected))=NaN;
        assert(isequaln(images(j).CData,expected));
        actual=actual+1;
    end
    close(f); assert(actual==2,'Near-bottom native panels missing');
    rows=struct([]); exact=true;
    for j=1:5
        file=fullfile(v,[stems{j} '.png']); info=imfinfo(file); assert(info.Width==1600 && info.Height==1000);
        f=openfig(fullfile(v,[stems{j} '.fig']),'new','invisible');
        dest=fullfile(a,[stems{j} '_roundtrip.png']); print(f,dest,'-dpng','-r100'); close(f);
        p=imread(file); q=imread(dest); assert(isequal(size(p),size(q)));
        maxdiff=max(abs(double(p(:))-double(q(:)))); changed=nnz(any(p~=q,3));
        exact=exact && maxdiff==0;
        rows=[rows struct('figure',string(stems{j}),'reopened',true,'width',info.Width,'height',info.Height, ...
            'pixel_max_difference',maxdiff,'changed_pixels',changed)]; %#ok<AGROW>
    end
    writetable(struct2table(rows),fullfile(a,'figure_roundtrip.csv'));
    receipt.status='FIGURE_DATA_COORDINATES_AND_EDITABILITY_PASS'; receipt.exit_status=0;
    receipt.required_checks=struct('full_map_actual_uniform_depths',true,'full_map_saved_tl_identical',true, ...
        'axis_original_data_identical',true,'near_bottom_original_data_identical',true,'png_dimensions',true,'figures_reopened',5);
    receipt.pixel_roundtrip_all_exact=exact;
    receipt.pixel_note='Strict pixel results recorded separately; scientific data and physical coordinate identity are required. No pixel threshold was relaxed.';
catch err
    if exist('f','var') && isgraphics(f), close(f); end
    receipt.status='FIGURE_AUDIT_FAILED'; receipt.exit_status=1; receipt.failed_checks=1;
    receipt.error=getReport(err,'extended','hyperlinks','off');
end
receipt.elapsed_s=toc(t); fid=fopen(fullfile(a,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,'PrettyPrint',true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt));
end
