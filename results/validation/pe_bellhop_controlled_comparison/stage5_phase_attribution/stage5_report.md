# Stage 5 phase discrepancy attribution

状态：**PASS**；solver calls: `0`。

代表点规则：first non-I height, maximum sampled height, maximum sampled curvature。分类严格使用 Goal 的 rho_shape 与 E_aligned/E_G 门限。

| A | K | region | E_G | E_aligned | reduction | phi0 | rho_shape | phase RMS | TL RMS | phase-gradient RMS/max | PE/BH vs stationary RMS | 2k eta vs angle RMS | classification | guards |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---:|
| 0.05 | 0.1 | II | 0.0210357376 | 0.0193849209 | 7.848% | -0.00818966942 | 0.99981227 | 0.0208690864 | 0.0228496268 | 0.00215310892/0.0055932975 | 0.0208684297/1.72948547e-06 | 0.0202833139 | spatial distortion | PASS |
| 0.2 | 0.1 | II | 0.100485227 | 0.0801428322 | 20.244% | -0.0609643815 | 0.996795793 | 0.100141661 | 0.0914359332 | 0.00810281755/0.0160911249 | 0.100139903/6.17857192e-06 | 0.0825710388 | spatial distortion | PASS |
| 0.02 | 0.47 | II | 0.0228160865 | 0.0210042807 | 7.941% | -0.00891676522 | 0.999779437 | 0.0129059305 | 0.163441183 | 0.00654944562/0.0261215013 | 0.0129018071/1.29065286e-05 | 0.00709721099 | spatial distortion | PASS |

Alignment removes less than 50% of E_G for every representative case, while rho_shape remains high. Under the frozen decision rule these cases are spatial-distortion dominated, not global/coherent-phase dominated.
The stationary/specular, phase-gradient, 2k eta, and angle-aware phase quantities are diagnostics only; no production PE screen or Bellhop physics is changed.
