function checks=assert_spectral_fast_kirchhoff_multifrequency_artifacts()
%ASSERT_SPECTRAL_FAST_KIRCHHOFF_MULTIFREQUENCY_ARTIFACTS Check saved evidence.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','spectral_fast_kirchhoff_multifrequency_consistency');
q=load(fullfile(out,'consistency.mat'),'s'); s=q.s;
[eta,sl,sec]=evaluate_pm_pilot_surface(s.profile,s.x_m);
assert(isequal(eta,s.eta_m) && isequal(sl,s.slope) && isequal(sec,s.second));
assert(numel(s.sparse)==5 && isequal([s.sparse.frequency_hz]',s.sparse_hz));
assert(all([s.sparse.rank_in]==s.ranks(1)) && all([s.sparse.rank_out]==s.ranks(2)));
assert(numel(s.anchors)==3 && numel(s.guards)==3);
for a=[s.sparse s.dense]
 assert(all(isfinite([a.full;a.fast])));
 m.whole=norm(a.fast-a.full)/norm(a.full);
 w=s.weights(s.mask); w=w/sum(w); aa=a.fast(s.mask); bb=a.full(s.mask);
 m.weighted=sqrt(sum(w.*abs(aa-bb).^2)/sum(w.*abs(bb).^2));
 assert(abs(m.whole-a.metrics.whole)<1e-14 && abs(m.weighted-a.metrics.weighted)<1e-14);
 k=2*pi*a.frequency_hz/1500;
 assert(abs(a.H_full_physical-a.H_full_reduced*exp(1i*k*103))<1e-14);
 % tau_ref-tau_direct and 6/c use different floating-point evaluation orders.
 assert(abs(a.H_full_dsp-a.H_full_reduced*exp(-1i*k*6))<1e-12);
 assert(abs(a.H_full_native_dsp-conj(a.H_full_physical)*exp(1i*k*97))<1e-14);
 bridge=apply_pe_channel_phase_reference_vertical(a.frequency_hz, ...
  struct('direct_f',a.H_direct_native_dsp,'reflect_fm',conj(a.H_full_reduced)),s.geometry,'direct_dsp');
 assert(abs(bridge.reflect_fm-a.H_full_native_dsp)<1e-12);
 if s.consistency_pass
  assert(abs(a.H_total_full-a.H_direct_native_dsp-a.H_full_native_dsp)<1e-14);
  assert(abs(a.H_total_fast-a.H_direct_native_dsp-a.H_fast_native_dsp)<1e-14);
 end
end
if numel(s.dense)==65
 assert(s.sparse_pass && isequal([s.dense.frequency_hz]',s.dense_hz));
 hf=complex(zeros(65,1)); ha=hf; hb=hf;
 for i=1:65
  hf(i)=s.dense(i).H_full_native_dsp; ha(i)=s.dense(i).H_full_dsp; hb(i)=s.dense(i).H_fast_dsp;
 end
 assert(norm(s.phase.cir_full-ifft(ha))<1e-14);
 assert(norm(s.phase.cir_fast-ifft(hb))<1e-14);
 assert(norm(s.native_phase.cir_full-ifft(hf))<1e-14);
 t=readtable(fullfile(out,'dense_native_physical_H.csv'));
 assert(norm(complex(t.Full_native_dsp_real,t.Full_native_dsp_imag)-hf)<1e-14);
 physical=exp(1i*s.native_phase.full_physical_unwrapped);
 for i=1:65, assert(abs(physical(i)-s.dense(i).H_full_physical/abs(s.dense(i).H_full_physical))<1e-11); end
 assert(max(abs(s.native_phase.full_dsp_group_delay_s+s.native_phase.full_physical_group_delay_s+97/1500))<1e-11);
end
checks=struct('surface_bit_identical',true,'complex_metrics_recomputed',true, ...
 'phase_reference_closure',true,'unified_rank',s.ranks,'sparse_count',numel(s.sparse), ...
 'dense_count',numel(s.dense),'scientific_status',s.status);
fid=fopen(fullfile(out,'artifact_assertions.json'),'w');
fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true)); fclose(fid);
disp(checks);
end
