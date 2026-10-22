# Optional Python plotting workflow

Read this file only for a Python plotting task or when a bundled Python audit helper is needed. MATLAB plotting does not require Python, uv, Matplotlib, or pypdf.

Generic Python exports default to PNG. PDF remains supported only when explicitly requested or included in the current task's delivery requirements. Vector-format examples and PDF font checks below apply only when that format is actually requested.

Commands and paths below are relative to the skill root. The pinned dependencies are a dated optional snapshot, not installation requirements for MATLAB.

### 4. Implement with scoped styles

Use Matplotlib's object-oriented API and temporary style contexts:

```python
import matplotlib.pyplot as plt

from style_presets import style_context

with style_context("default", palette_name="okabe_ito_on_white"):
    fig, ax = plt.subplots(
        figsize=(89 / 25.4, 60 / 25.4),
        layout="constrained",
    )
    ax.plot(x, y, marker="o", label="Observed")
    ax.set(xlabel="Time (hours)", ylabel="Response (unit)")
    ax.legend()
```

`layout="constrained"` supports colorbars, nested GridSpec, subfigures, and `subplot_mosaic`. Do not call `tight_layout()` afterward; it disables constrained layout.

For exact physical dimensions, do not use `bbox_inches="tight"` unless the changed page size is intentional.

#### Color normalization

```python
import matplotlib as mpl

norm = mpl.colors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=5)
cmap = mpl.colormaps["RdBu_r"].with_extremes(bad="#777777")
image = ax.imshow(values, norm=norm, cmap=cmap, interpolation="nearest")
fig.colorbar(image, ax=ax, label="Change (unit)")
```

Use `LogNorm`, `CenteredNorm`, `SymLogNorm`, `BoundaryNorm`, or `TwoSlopeNorm` only when its mapping matches the scientific meaning.

#### Seaborn

Seaborn 0.13.2 uses the current `errorbar` API:

```python
sns.lineplot(
    data=frame,
    x="time",
    y="response",
    hue="treatment",
    style="treatment",
    markers=True,
    errorbar=("ci", 95),
    n_boot=5000,
    seed=20260723,
    ax=ax,
)
```

Axes-level functions fit custom Matplotlib layouts; figure-level functions create their own figures/facets. Do not customize Seaborn's internal artist lists as if they were stable API.

#### Plotly

- Use `write_html()` for interaction and `write_image()`/`plotly.io.write_images()` for static output.
- Kaleido 1.3.0 requires Chrome/Chromium; it no longer bundles Chrome.
- Current static formats: PNG, JPEG, WebP, SVG, PDF. EPS is Kaleido v0-only.
- Do not pass deprecated `engine=` or use Orca/`plotly.io.kaleido.scope`.
- `width`, `height`, and `scale` control pixels; `scale=3` is not inherently “300 DPI.”
- WebGL traces embed raster content in PDF/SVG.
- Fully offline exports need local external assets when a figure references MathJax/topojson/tiles.

### 5. Export explicitly and record provenance

```python
from figure_export import export_figure

report = export_figure(
    fig,
    "outputs/figure1",
    formats=["png"],
    dpi=600,
    bbox_inches=None,  # preserve figure page dimensions
    provenance={
        "raw_data": "data/source.csv",
        "transformations": ["predeclared QC filter", "group mean"],
        "uncertainty": "95% bootstrap CI; seed 20260723",
        "missing_data": "retained as gaps",
    },
    write_manifest=True,
)
```

The exporter refuses implicit overwrite, writes atomically, keeps vector DPI for embedded rasters, uses TIFF LZW, and can use PDF/PS Type 42 fonts. It does not validate scientific content or publisher acceptance.

For editable fonts:

- PDF/PS Type 42 embeds TrueType fonts.
- `svg.fonttype="none"` keeps text editable/searchable but does not embed fonts; appearance depends on installed fonts.
- `svg.fonttype="path"` preserves glyph appearance as paths but loses editable/searchable text.

Use an opaque explicit background unless transparency is required; blending against another background changes apparent contrast.

### 6. Inspect, compare, and review

1. Inspect file metadata.
2. Audit palette contrast/grayscale separation.
3. Compare against a dated publisher snapshot.
4. View at final size in the manuscript/web context.
5. Manually review fonts, embedded rasters, clipping, legends, scale bars, image integrity, caption, alt text, and source data.
6. Re-check the live target-journal page immediately before upload.

## Pinned snapshot

The examples and smoke tests use direct package pins current on 2026-07-23:

```bash
uv run --isolated --no-project --python 3.13 \
  --with "matplotlib==3.11.1" \
  --with "seaborn==0.13.2" \
  --with "plotly==6.9.0" \
  --with "kaleido==1.3.0" \
  --with "pillow==12.3.0" \
  --with "pypdf==6.14.2" \
  python your_figure.py
```

This is a dated direct-dependency snapshot, not a transitive lock. Use the project's uv lock for exact replay; this skill intentionally ships no dependency lock.

## Bundled CLIs

All helpers are deterministic, network-free, bounded, reject symlink inputs/destinations where relevant, and refuse overwrite unless `--force` is explicit.

### Inspect raster/vector metadata

```bash
uv run --isolated --no-project --python 3.13 \
  --with "pillow==12.3.0" \
  python scripts/image_metadata.py figure.tiff \
  --format tiff --mode RGB --min-dpi 300 --target-width-mm 85 \
  --alpha-policy forbid
```

Supports raster images (Pillow), SVG, PDF (pypdf), and EPS/PS. Reports dimensions, DPI/effective DPI, mode, alpha, ICC presence, compression, page size, and conservative first-page PDF font resources. It does not inspect every embedded raster in a vector container.

### Audit palette contrast and grayscale

```bash
uv run --isolated --no-project --python 3.13 \
  python scripts/palette_audit.py \
  --palette okabe_ito_on_white \
  --background FFFFFF \
  --role graphical
```

Reports exact WCAG sRGB contrast plus pairwise CIE L* grayscale screening. The grayscale threshold is a heuristic, not a standard.

### Plan/screen publisher export

```bash
uv run --isolated --no-project --python 3.13 \
  python scripts/export_plan.py \
  --publisher nature \
  --figure-type combination \
  --width single \
  --phase final
```

Add `--input figure.pdf` to screen machine-readable properties. Profiles are official-source snapshots accessed 2026-07-23, not automatic compliance rules.

### Preview styles

```bash
uv run --isolated --no-project --python 3.13 \
  --with "matplotlib==3.11.1" \
  python scripts/style_preview.py \
  --output outputs/style-preview \
  --style default \
  --palette okabe_ito_on_white \
  --formats png,svg
```

### Inspect/write styles and smoke-test export

```bash
uv run --isolated --no-project --python 3.13 \
  python scripts/style_presets.py --list
uv run --isolated --no-project --python 3.13 \
  python scripts/style_presets.py --show nature
uv run --isolated --no-project --python 3.13 \
  --with "matplotlib==3.11.1" \
  python scripts/figure_export.py --demo outputs/export-smoke --manifest
```

## Assets

- `assets/publication.mplstyle`: general print starting point.
- `assets/nature.mplstyle`: dated flagship Nature visual starting point, not a compliance preset.
- `assets/presentation.mplstyle`: larger projected-display style.
- `assets/color_palettes.py`: importable Okabe-Ito and Paul Tol values with metadata.
- `assets/publisher_profiles.json`: dated, machine-readable planning snapshots.

Matplotlib style files omit `#` in hex colors because `#` begins comments in `.mplstyle` parsing.

