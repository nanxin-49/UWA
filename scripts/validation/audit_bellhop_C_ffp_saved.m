function receipt=audit_bellhop_C_ffp_saved(project_root,comparison_dir,figure_dir)
% Independent frozen-source, saved-array, CSV and raw-field SHA audit.
start=tic; checks=table(strings(0,1),false(0,1),'VariableNames',{'check','pass'}); s=load(fullfile(comparison_dir,'comparison.mat')); meta=s.meta;
a=load(meta.C_file,'total','envelope'); b=load(meta.FFP_file,'output'); add('C_pressure_array_exact_frozen',isequaln(s.p_bellhop,a.total)); add('FFP_pressure_array_exact_frozen',isequaln(s.p_ffp,b.output.pressure_1m)); add('C_envelope_exact_frozen',isequaln(s.envelope,a.envelope)); clear a b
add('C_frozen_SHA_after',strcmpi(black_file_hash(meta.C_file),meta.C_sha256)); add('FFP_frozen_SHA_after',strcmpi(black_file_hash(meta.FFP_file),meta.FFP_sha256));
base=fullfile(project_root,'results','validation','bellhop_rotated_fig12b'); fields={fullfile(base,'20261115_A_resume_01','assembly','result.mat'),fullfile(base,'20261115_B_run_02','assembly','result.mat')}; expected={'e94e79fc6e235ecf36db2b6bdbd88bc0911f98d5114fa0c364b32b13f89cc6b4','752062dfd553fc05abb1b568e2e209862810df5ea84a8230c3735ba0a2d3ea62'};
for k=1:2, add(sprintf('AB_field_SHA_%d',k),strcmpi(black_file_hash(fields{k}),expected{k})); end
t=readtable(fullfile(comparison_dir,'common_point_metrics.csv'),'Delimiter',',','TextType','string','VariableNamingRule','preserve'); [z,x]=ndgrid(s.z,s.x); r=hypot(x,z-4950); valid=~s.source_mask; regions={valid,valid&z<=4900,valid&z>=4900&r>=1,valid&x==0&r>=1,valid&r<1,valid&abs(x)<=1};
for k=1:numel(regions)
 use=regions{k}; d=s.p_bellhop(use)-s.p_ffp(use); rel=norm(d)/norm(s.p_ffp(use)); mx=max(abs(d)./s.envelope(use)); support=use&abs(s.p_bellhop)>=.05*s.envelope&abs(s.p_ffp)>=.05*s.envelope; ph=max(abs(angle(s.p_bellhop(support).*conj(s.p_ffp(support)))));
 add(sprintf('region%d_complexL2_CSV',k),abs(rel-t.complex_relative_L2_FFP(k))<1e-14); add(sprintf('region%d_max_over_envelope_CSV',k),abs(mx-t.max_complex_over_C_envelope(k))<1e-14); add(sprintf('region%d_phase_CSV',k),abs(ph-t.max_effective_phase_rad(k))<1e-14); add(sprintf('region%d_support_CSV',k),nnz(support)==t.phase_support(k));
end
add('only_true_source_mask',nnz(s.source_mask)==1&&isequal(isnan(s.p_bellhop),s.source_mask)&&isequal(isnan(s.p_ffp),s.source_mask));
figures=dir(fullfile(figure_dir,'F*.fig')); add('all_six_native_FIG',numel(figures)==6); for k=1:numel(figures), [~,stem]=fileparts(figures(k).name); info=imfinfo(fullfile(figure_dir,[stem '.png'])); add([stem '_PNG'],info.Width>1000&&info.Height>500); end
receipt=struct('status','FROZEN_SOURCE_COMMON_POINT_INDEPENDENT_AUDIT','checks_count',height(checks),'failed_checks',nnz(~checks.pass),'elapsed_s',toc(start),'new_native',0,'new_FFP',0,'ABC_recomputed',0,'numerical_acceptance','Artifact/data consistency only; no new cross-method physical threshold'); writetable(checks,fullfile(comparison_dir,'independent_checks.csv')); fid=fopen(fullfile(comparison_dir,'checkpoint_independent.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); fprintf('C_FFP_INDEPENDENT checks%d failed%d elapsed%.3fs\n',height(checks),receipt.failed_checks,receipt.elapsed_s);
 function add(name,pass), checks(end+1,:)={string(name),logical(pass)}; end
end
