# Spectral Fast FK multi-frequency consistency

Status: **PASS**. Date: 2026-10-08. MATLAB 25.2.0.3312555 (R2025b) Update 6.

本轮 multi-frequency consistency **PASS**：固定 seed260001/C2 海面在4–8 kHz的65个频点上，统一3/12 rank满足FastFK–FullFK复误差门槛，三个BIE物理锚点、数值护栏、相位/群时延/CIR一致性及加速门均通过。冻结3/8在6–8 kHz失败，3/10在8 kHz仍略超标，因此采用整个频带统一3/12，未做逐频tuning。

相位桥接不能只补载波。原生Full/BIE/Fast物理相量到MATLAB IFFT传递函数采用一次明确的全包络Fourier符号转换：conj(P)*exp(+iωτdir0)，三条路径共享同一个直达时刻参考。验证层经该数学转换后再使用现有direct-DSP helper；这不是逐频相位拟合，也没有改动reference或生产API。未转换的carrier-only导入仍保留为失败的独立控制。后续固定PM宽带及小规模realization数据可以沿这个明确记录的表示继续；本轮不把原生FK reduced字段直接当作公共PE reduced字段或宣称生产集成已完成。

Canonical seed260001; 12 immutable Fourier modes; actual Kmax 0.471238898 rad/m; coefficient SHA256 `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`. Surface x=[-60,60] m, N4097; C2 plateau42/support50 m, common height bound 0.5 m. All frequencies use bit-identical eta, analytic slopes and coefficients, no height scaling, complex scalar fitting or phase alignment. Source: frozen984-point Gaussian sigma=.3 m, period192.1875 m, x0=-96.09375 m; all saved bins retained, frequency-dependent exact kz. Uniform c1500, Tx100/Rx3, 2-D line-source validation chain; production PE is unchanged.

Sparse axis4/5/6/7/8 kHz; BIE anchors4/6/8 kHz; selected common rank3/12, padding2. Dense follow-up65 samples, delta-f62.5 Hz, window16 ms, bin spacing16/65=.246154 ms; bandwidth4 kHz gives nominal resolution about.25 ms. Sparse delta-f1 kHz has only1 ms window and aliases nominal4 ms relative delay: no sparse CIR or physical group-delay certification. Dense execution requires sparse algorithm, physics, numerical and cold speed gates.

The16 ms window contains the nominal4 ms reflection excess delay, but does not certify all scattering tails as unaliased. A point-source geometric path through x=60 m at the flat integration edge permits excess delay(sqrt(100^2+60^2)+sqrt(3^2+60^2)-97)/1500=53.129 ms. Such contributions may be weak/nonstationary; their absence has not been proved. CIR remains a finite-band periodic diagnostic, not a support-converged final response. Samples are IFFT-ordered as frequency offsets from4 kHz (no ifftshift); a different baseband carrier changes complex modulation, not the compared magnitude/timing.

## Phase reference

Physical exp(-i omega t) Full/Fast outputs carry exp(+ik*103). Native reduced reflection is P*exp(-ik*103); it is not automatically treated as the public PE reduced field. The common-origin MATLAB IFFT transfer is **H_native_dsp=conj(P)*exp(+ik*97)**. Equivalently, supply conj(native reduced) to the unchanged `apply_pe_channel_phase_reference_vertical` helper for its direct_dsp target, which multiplies by exp(-ik*6). The direct source uses the same saved Gaussian coefficients propagated97 m; its DSP envelope is conj(Pdirect*exp(-ik*97)). This is one deterministic Fourier-sign conversion, not a fitted complex scalar, per-frequency phase zero or peak alignment. Native physical fields remain untouched and are not reconstructed by calling the helper absolute_physical target on a DSP-adapted envelope. Total is formed only after reflected-only consistency gates. Public PE/interface semantics are unchanged.

Raw physical phase is saved; blind unwrap is aliased even at62.5 Hz because nominal absolute flight time68.667 ms exceeds16 ms. Carrier-restored physical unwrap is unwrap(native reduced phase)+2*pi*f*103/1500. Its absolute branch is only defined modulo2*pi. Negative -dphi/(2*pi*df) for a positive physical carrier is expected; the IFFT-ready response has the opposite Fourier sign.

### Independent native-physical to IFFT bridge control

The unadapted candidate exp(-ik*6)*(P*exp(-ik*103)) is retained only as a carrier-only import diagnostic (H_full_dsp/H_fast_dsp in MAT and dense_H.csv). It is not the primary physical transfer. A fixed-kx4 rad/m plane wave P=-exp(i*kz*103) gives native relative delay=(103*k/kz-97)/1500; the unadapted import gives6/1500-103*(k/kz-1)/1500. They agree at kx0 and differ for kx nonzero: maximum sampled delay difference 4.023489 ms. After the whole-envelope conversion, helper output matches the independently constructed native physical IFFT transfer at relative complex L2 2.067614e-13; bridge pass=1. This validates the external native FK import, not a production API migration. Primary tables are dense_native_physical_H.csv, phase_group_delay.csv and reflected_cir.csv; phase_reference_meta documents the source/target representation and single time origin.

## Per-frequency algorithm fidelity and descriptors

|Hz|weighted L2|whole L2|magnitude|phase RMS rad|effective in/out|Kmax/k|curvature RMS/k|illuminated curvature/k|out theta99 deg|outside incident99|
|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
|4000|0.00023066996|0.00028748516|8.688041e-05|0.0002397295|3/7|0.028125|0.001192235|0.001477154|20.86474|0.005758169|
|5000|0.00011077283|0.00014125447|3.282947e-05|0.0001772396|3/8|0.0225|0.0009537876|0.001205882|19.71962|0.05691989|
|6000|5.879272e-05|9.2574575e-05|1.597215e-05|0.0001593824|3/9|0.01875|0.000794823|0.001023427|19.087|0.1835906|
|7000|3.728881e-05|6.1131037e-05|1.042815e-05|0.0004496241|3/10|0.01607143|0.0006812769|0.0008934055|18.74459|0.3045354|
|8000|4.0199509e-05|5.2331609e-05|2.423514e-05|0.001002219|3/10|0.0140625|0.0005961173|0.0007956824|18.3939|0.4609813|

Fixed canonical Stage0 M99 weights/mask are held across frequencies; whole-grid error is separate. All-native complex fields are retained. Curvature is eta_second/(1+slope^2)^(3/2), illumination weights are incident reference-plane energy at each frequency. Spectrum is a finite Hann-aperture diagnostic, not far-field ocean flux; no new universal validity law.

## Physical anchors

|Hz|Full--BIE weighted|whole line diagnostic|phase RMS|boundary residual|nodes|wall s|matrix MiB|
|---:|---:|---:|---:|---:|---:|---:|---:|
|4000|0.00080009587|0.009965395|0.0007857323|7.690677e-08|2560|150.72|100|
|6000|0.00053037876|0.0005884464|0.0005104319|4.11751e-08|3840|284.86|225|
|8000|0.00040826892|0.0004019282|0.0003829177|5.062065e-08|5120|481.05|400|

BIE uses L60, ppw8/order8/self96, h=-2, plateau.85, the same analytical C2 surface and incident angular spectrum. Physical gate is weighted error<1e-3 and offgrid residual<1e-4, restricted to the frozen M99 core; edge errors are diagnostic. Numerical BIE refinement is required if a residual or near-threshold ambiguity appears.

## Numerical guards

|Hz|padding2/4|height nodes33/65|Full4097/8193 weighted|whole|Fast grid weighted|eta closure|
|---:|---:|---:|---:|---:|---:|---:|
|4000|4.006815e-16|2.193067e-10|6.560588e-06|1.629997e-05|6.562179e-06|0|
|6000|4.790245e-16|3.332706e-08|2.353252e-07|5.774076e-07|2.388994e-07|0|
|8000|4.722712e-16|1.270575e-06|1.090013e-07|1.06581e-07|1.336132e-06|0|

Execution note: the8 kHz BIE exceeded the MCP300-second response wait; the original shared-session computation continued and saved its checkpoint after about481 seconds. A subsequent harmless MCP probe confirmed recovery. No duplicate BIE run or shell MATLAB launch was used.


## CPU and memory

|Hz|Full setup/apply wall s|Fast setup/apply wall s|warm/cold speed|plan MiB|setup/apply arrays MiB|Full block estimate MiB|
|---:|---|---|---|---:|---|---:|
|4000|0.12311 / 2.2438|0.16807 / 0.036903|60.8 / 11.55|8.3195|29.037 / 0.87521|104.03|
|5000|0.098437 / 2.1655|0.15699 / 0.030656|70.64 / 12.06|8.3195|29.037 / 0.87521|104.03|
|6000|0.21165 / 2.6054|0.15702 / 0.030045|86.71 / 15.06|8.3195|29.037 / 0.87521|104.03|
|7000|0.15708 / 2.5853|0.14672 / 0.029851|86.61 / 15.53|8.3195|29.037 / 0.87521|104.03|
|8000|0.14202 / 2.4821|0.14968 / 0.031539|78.7 / 14.48|8.3195|29.037 / 0.87521|104.03|

Fast apply median of3 after warmup (selected sparse rank remeasured with timing mode); Full one apply per frequency, timing is hardware/session-dependent. CPU columns are in metrics CSV. Memory is named arrays/analytical block estimate, not measured process peak. Geometry eta/slopes/normal/quadrature and Gaussian coefficients are per-realization reusable; kz, both SVD bases, height carrier, Hankel lag kernels and propagation phases are per-frequency. Existing frozen implementation rebuilds them; no claimed cross-frequency setup acceleration. Geometry reuse is already represented by fixed surface state; its negligible load/setup time is not included in cold speed.

Measured one-time fixed Fourier geometry evaluation/normal/quadrature setup: 0.0032177 wall s, 0.015625 CPU s, named geometry arrays 0.1875458 MiB. Saved-source/input file IO is outside this kernel benchmark.

## Dense phase and CIR diagnostic

All65-point algorithm gate: 1. Max native relative phase 0.00017439397 rad; max DSP group-delay difference 7.0193722e-08 s; CIR relative complex L2 7.1717713e-05; max adjacent reduced phase 0.1591546 rad. Sparse/dense boundary is governed by fixed gates, never total-channel cancellation. IFFT uses the rectangular65-sample native-physical direct-DSP reflected-only response with no taper, zero padding or peak alignment; diagnostic only, not a converged final CIR.

Across65 samples: max weighted/whole-grid complex errors 0.000230669959 / 0.000287485157; minimum cold/warm speedups 9.832663 / 41.93941. Native direct-DSP group delay range 3.5947172--3.59507788 ms; sampled CIR maximum at 3.69230769 ms. The grid bin size explains why the sampled peak differs from the approximately3.595 ms phase delay; no peak alignment is performed.

The unadapted carrier-only diagnostic also preserves Full/Fast fidelity, but that alone does not establish its native physical delay interpretation. All primary H/CIR results use the validated native-to-DSP bridge and one shared direct time origin.


## Frozen-rank failure and common-rank decision

|Hz|3/8 weighted error|3/8 whole error|
|---:|---:|---:|
|4000|0.00024763757|0.00029660942|
|5000|0.00052823153|0.00053293787|
|6000|0.0021891266|0.0021813652|
|7000|0.0074891234|0.0074446234|
|8000|0.022496677|0.022316181|

The frozen3/8 fails from6 kHz. Common3/10 narrowly fails at8 kHz; incoming4/5 at outgoing10 does not resolve it. The next outgoing ladder3/12 passes all sparse frequencies. All trials are saved in uniform_rank_trials.csv; no frequency-specific tuning. Chosen rank3/12 applies to the complete sparse/dense band.

## Decision and evidence boundary

Algorithm=1; physics=1; numerical=1; speed=1; sparse=1; multi-frequency PASS=1. Full/Fast frequency consistency=1; native-to-project phase bridge=1. Allows the next separately scoped fixed-PM broadband and small realization data stage under these tested conditions only.

Verification: five independent tests cover nominal4 ms delay, sparse aliasing, fixed surface nesting, nonzero-kx native/DSP conversion, and the actual phase assembler with a known4.3 ms response. Artifact assertions recompute complex errors from raw fields, compare CSV imaginary signs to scalar native fields, verify exp(i*restored physical phase), check the physical/DSP derivative sign identity, and confirm total-component closure. Code Analyzer is clean for the four new MATLAB files. Complex frequency collection uses vertcat (no implicit conjugate transpose). Reference-phase algebra checks use1e-12 to cover alternative evaluation orders of large carrier arguments; mathematical/scientific fidelity gates remain1e-3 and the delay-difference gate1 microsecond. No production or communication regression run was needed because this task only adds validation workflows and explicitly excludes communication work.

Artifacts: `results/validation/spectral_fast_kirchhoff_multifrequency_consistency/consistency.mat`, CSV tables, raw complex fields, fixed surface, input/reference SHA256 provenance and figures. Reference FullFK/BIE/FastFK implementations and public paramsV/output remain unchanged.

![Sparse audit](../results/validation/spectral_fast_kirchhoff_multifrequency_consistency/figures/sparse_audit.png)

![Dense phase and CIR](../results/validation/spectral_fast_kirchhoff_multifrequency_consistency/figures/dense_phase_cir.png)

![Native physical IFFT and phase bridge control](../results/validation/spectral_fast_kirchhoff_multifrequency_consistency/figures/native_ifft_reference_control.png)
