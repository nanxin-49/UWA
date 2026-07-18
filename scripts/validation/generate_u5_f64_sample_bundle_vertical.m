function meta=generate_u5_f64_sample_bundle_vertical(n_samples,seed)
%GENERATE_U5_F64_SAMPLE_BUNDLE_VERTICAL Save H and physical CIR without PE.
if nargin<1, n_samples=10000; end
if nargin<2, seed=2700001; end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
folder=fullfile(root,'results','validation','u5_conditional_channel_f64');
S=load(fullfile(folder,'u5_conditional_channel_model_f64_full.mat'),'model'); model=S.model;
timer=tic; [draw,sample_meta]=sample_conditional_channel_vertical(model,n_samples,seed, ...
    struct('path','auto','rank','selected')); sample_elapsed=toc(timer);
timer=tic; cir=build_physical_cir_vertical(draw.H_total_f,model.frequency_axis, ...
    model.reference_delay_s,'none',1); cir_elapsed=toc(timer);
samples=draw; samples.h_physical_tau=cir.h_physical_tau;
samples.delay_axis=cir.delay_axis_s; samples.frequency_axis=model.frequency_axis;
samples.reference_delay_s=model.reference_delay_s;
samples.physical_delay_resolution_s=cir.physical_delay_resolution_s;
samples.maximum_unambiguous_delay_s=cir.maximum_unambiguous_delay_s;
samples.zero_padding_note=cir.zero_padding_note;
samples.model_kind=model.kind;
meta=struct('n_samples',n_samples,'seed',seed,'sample_elapsed_s',sample_elapsed, ...
    'cir_elapsed_s',cir_elapsed,'total_elapsed_s',sample_elapsed+cir_elapsed, ...
    'per_channel_s',(sample_elapsed+cir_elapsed)/n_samples, ...
    'sampling_meta',sample_meta,'pe_calls',0);
file=fullfile(folder,sprintf('u5_f64_generated_samples_%d.mat',n_samples));
timer=tic; save(file,'samples','meta','-v7.3'); save_elapsed=toc(timer);
info=dir(file); meta.file_bytes=info.bytes; meta.save_elapsed_s=save_elapsed;
meta.total_with_save_s=meta.total_elapsed_s+save_elapsed;
save(fullfile(folder,sprintf('u5_f64_generated_samples_%d_meta.mat',n_samples)),'meta');
fprintf('Generated %d H/CIR samples in %.6f s, save %.6f s, total %.6f s, PE calls=%d, file %.3f MiB\n', ...
    n_samples,meta.total_elapsed_s,meta.save_elapsed_s,meta.total_with_save_s,meta.pe_calls,info.bytes/2^20);
end
