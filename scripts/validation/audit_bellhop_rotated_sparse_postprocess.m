function receipt = audit_bellhop_rotated_sparse_postprocess(comparison_dir,figure_dir,output_dir,options)
%AUDIT_BELLHOP_ROTATED_SPARSE_POSTPROCESS Saved arithmetic/CSV/FIG/PNG replay.
% This audits the sparse diagnostic, not complete-grid or physical validity.
arguments
    comparison_dir (1,:) char
    figure_dir (1,:) char
    output_dir (1,:) char
    options.reference_subset_file (1,:) char=''
    options.subset_manifest_file (1,:) char=''
end
assert(~isfolder(output_dir)); [~,mem]=memory;
needed=2^30; if ~isempty(options.reference_subset_file), needed=512*2^20; end
assert(mem.PhysicalMemory.Available>=needed,'Independent saved-data audit RAM gate.');
mkdir(output_dir); t=tic; diary(fullfile(output_dir,'diary.log')); cleaner=onCleanup(@()diary('off'));
d=load(fullfile(comparison_dir,'comparison.mat')); assert(strcmp(d.receipt.scope,'pilot'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
assert(strcmpi(black_file_hash(d.receipt.field_file),d.receipt.field_sha256));
assert(strcmpi(black_file_hash(d.receipt.ffp_file),d.receipt.ffp_sha256));
source=load(d.receipt.field_file,'total','envelope');
add('raw_field_binding',max(abs(source.total-d.p_bellhop),[],'all','omitnan'),0);
add('raw_envelope_binding',max(abs(source.envelope-d.envelope),[],'all','omitnan'),0);
reference_mode='FULL_REFERENCE_READ';
if isempty(options.reference_subset_file)
    f=load(d.receipt.ffp_file,'output');
    [yes,iz]=ismembertol(d.z,f.output.depth_m(:),1e-9,'DataScale',1); assert(all(yes));
    [yes,ix]=ismembertol(d.x,f.output.range_m(:).',1e-9,'DataScale',1); assert(all(yes));
    ref=f.output.pressure_1m(iz,ix);
else
    assert(isfile(options.subset_manifest_file),'Validated subset requires frozen artifact manifest.');
    manifest=jsondecode(fileread(options.subset_manifest_file));
    entry=find(strcmpi({manifest.path},options.reference_subset_file)); assert(isscalar(entry));
    assert(strcmpi(black_file_hash(options.reference_subset_file),manifest(entry).sha256),'Independent subset SHA mismatch.');
    f=load(options.reference_subset_file,'reference','comparison','plan');
    assert(strcmpi(f.comparison.ffp_sha256,d.receipt.ffp_sha256));
    assert(isequal(f.plan.depths_m(:),d.z) && isequal(f.plan.x_m(:).',d.x));
    ref=f.reference; reference_mode='FROZEN_INDEPENDENT_SUBSET_REPLAY_NO_LARGE_REFERENCE_LOAD';
end
add('independent_FFP_subset',max(abs(ref-d.p_ffp),[],'all','omitnan'),0);
add('source_mask_both_fields',double(~isequal(isnan(d.p_bellhop),isnan(ref),d.source_mask)),0);
valid=~d.source_mask; delta=d.p_bellhop-ref;
raw_error=abs(delta)./source.envelope;
add('pointwise_envelope_error',max(abs(raw_error(valid)-d.envelope_error(valid))),1e-14);
phase_use=valid & min(abs(d.p_bellhop),abs(ref))>=source.envelope/20;
ratio=d.p_bellhop(phase_use)./ref(phase_use); ph=atan2(imag(ratio),real(ratio));
add('phase_by_independent_ratio',max(abs(ph-d.phase_difference(phase_use))),1e-12);
add('phase_mask_binding',double(~isequal(phase_use,d.phase_valid)),0);
ze=[d.z(1);(d.z(1:end-1)+d.z(2:end))/2;d.z(end)];
xe=[d.x(1);(d.x(1:end-1).'+d.x(2:end).')/2;d.x(end)];
w=diff(ze)*diff(xe).'; [X,Z]=meshgrid(d.x,d.z); R=hypot(X,Z-4950);
masks={valid,valid & Z<=4900,valid & Z>=4900 & R>=1,valid & X==0 & R>=1,valid & R<1};
for k=1:height(d.metrics)
    use=masks{k}; name=char(d.metrics.region(k));
    l2=sqrt(sum(real(delta(use)).^2+imag(delta(use)).^2)/sum(abs(ref(use)).^2));
    weighted=sqrt(sum(w(use).*abs(delta(use)).^2)/sum(w(use).*abs(ref(use)).^2));
    add([name '_point_count'],abs(nnz(use)-d.metrics.points(k)),0);
    add([name '_complex_L2'],abs(l2-d.metrics.complex_relative_L2(k)),1e-12);
    add([name '_weighted_complex_L2'],abs(weighted-d.metrics.cell_weighted_complex_relative(k)),1e-12);
    add([name '_maximum_error_E'],abs(max(raw_error(use))-d.metrics.max_envelope_error(k)),1e-12);
end
csv=readtable(fullfile(comparison_dir,'reference_metrics.csv'),'TextType','string');
add('CSV_region_names',double(~isequal(csv.region,d.metrics.region)),0);
add('CSV_numeric_binding',max(abs(csv{:,2:8}-d.metrics{:,2:8}),[],'all','omitnan'),1e-12);
it=find(d.z==0,1); top=max(abs(d.p_bellhop(it,:))./source.envelope(it,:));
add('top_residual_recomputed',abs(top-d.receipt.boundary.top_bellhop_max_envelope),1e-14);
add('top_failure_preserved',double(d.receipt.boundary.top_sampled_criterion_met),0);
add('sparse_resolution_not_promoted',double(d.receipt.resolution.point01m_requirement_met),0);
add('bottom_derivative_not_invented',double(~strcmp(d.receipt.boundary.bottom_status,'NOT_EVALUABLE_SPARSE_POINTS')),0);
add('raw_visual_copy_SHA',double(~strcmpi(black_file_hash(fullfile(figure_dir,'raw_comparison.mat')),black_file_hash(fullfile(comparison_dir,'comparison.mat')))),0);
clear f ref source
stems={'F01_sparse_receivers_and_limits','F02_sparse_axis_complex_comparison'};
for k=1:2
    fig=openfig(fullfile(figure_dir,[stems{k} '.fig']),'new','invisible'); close_fig=onCleanup(@()close(fig));
    if k==1
        obj=findobj(fig,'Tag','sampled_envelope_error'); assert(isscalar(obj));
        use=~d.source_mask;
        add('F01_actual_X',max(abs(obj.XData(:)-X(use))),0);
        add('F01_actual_depth',max(abs(obj.YData(:)-Z(use))),0);
        add('F01_actual_error_color',max(abs(obj.CData(:)-raw_error(use))),1e-14);
        ax=ancestor(obj,'axes'); add('F01_color_limits',max(abs(ax.CLim-[0 .04])),0);
        add('F01_depth_positive_down',double(~strcmp(ax.YDir,'reverse')),0);
    else
        tags={'axis_abs_bellhop','axis_abs_ffp','axis_phase_difference','axis_envelope_error'};
        col=find(d.x==0,1);
        expected={abs(d.p_bellhop(:,col)),abs(d.p_ffp(:,col)),d.phase_difference(:,col),d.envelope_error(:,col)};
        for n=1:numel(tags)
            obj=findobj(fig,'Tag',tags{n}); assert(isscalar(obj));
            add([tags{n} '_x'],max(abs(obj.XData(:)-d.z)),0);
            add([tags{n} '_y'],double(~isequaln(obj.YData(:),expected{n})),0);
            add([tags{n} '_no_unsampled_connections'],double(~strcmp(obj.LineStyle,'none')),0);
        end
    end
    file=fullfile(output_dir,[stems{k} '_replay.png']); exportgraphics(fig,file,'Resolution',120);
    a=imread(fullfile(figure_dir,[stems{k} '.png'])); b=imread(file);
    add([stems{k} '_PNG_dimensions'],double(~isequal(size(a),size(b))),0);
    if isequal(size(a),size(b)), add([stems{k} '_PNG_pixel_max'],max(abs(double(a(:))-double(b(:)))),0); end
    clear close_fig
end
html=fileread(fullfile(figure_dir,'index.html'));
links=regexp(html,'(?:href|src)="([^"#]+)"','tokens'); missing=0;
for k=1:numel(links), missing=missing+~isfile(fullfile(figure_dir,links{k}{1})); end
add('offline_HTML_local_links',missing,0);
add('no_unrequested_PDF',numel(dir(fullfile(figure_dir,'*.pdf'))),0);
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status','COMPLETE_SPARSE_POSTPROCESS_AUDIT_ONLY','failed_checks',nnz(~checks.pass), ...
    'checks_count',height(checks),'elapsed_s',toc(t),'result_dir',output_dir, ...
    'reference_mode',reference_mode,'minimum_available_bytes',needed, ...
    'scope','Raw FFP/field/metric/CSV, editable FIG and PNG replay; no complete-field/physics qualification', ...
    'comparison_sha256',black_file_hash(fullfile(comparison_dir,'comparison.mat')));
save(fullfile(output_dir,'result.mat'),'receipt','checks');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end