# PE--Bellhop PM Stage 1B dimensionality sensitivity

状态：**PASS_WITH_LIMITS**

固定条件：4 kHz、uniform c=1500 m/s、seed=260001、W=192.1875 m、nx=984、PE step=0.05 m、y-invariant eta(x,y)=eta(x)、on-axis receiver。1-transverse PE 与 production 2-transverse PE 的 reflected-only rough/flat ratio 仅用于量化 dimensionality sensitivity。

## Results

| ny | G_1T | G_2T | delta TL (dB) | delta phase (rad) | complex relative error |
|---:|---:|---:|---:|---:|---:|
| 256 | -0.0868692441183+0.988068678385i | -0.0974201037782+0.969019755976i | -0.15885546 | 0.012505279 | 0.021953988 |
| 512 | -0.0868692441183+0.988068678385i | -0.0974244073654+0.969018832201i | -0.15885982 | 0.01250977 | 0.021956905 |

2T ny endpoint sensitivity: -4.3580248e-06 dB, 4.4916211e-06 rad. The 2T-vs-1T difference is dimensional/source mapping diagnostic, not a Bellhop or PE implementation failure.

## Checks

- profile_provenance: PASS
- one_transverse_finite: PASS
- two_transverse_finite: PASS
- ny_sensitivity_finite: PASS
- all: PASS

## Interpretation

The surface samples are copied identically in y, and the production PE core is unchanged. This report quantifies the 1T-to-2T bridge contribution that must be reported alongside the PE/Bellhop Stage 1A model discrepancy.
