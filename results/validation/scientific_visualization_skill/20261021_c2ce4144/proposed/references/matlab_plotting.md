# MATLAB plotting and editable figures

Use MATLAB for this repository's plots and the existing MATLAB MCP session for execution. Use the connected release's supported APIs; do not assume the release mentioned in another skill is the installed release.

## Output contract

For every MATLAB figure retained or delivered, default to `name.fig` and `name.png`. Use `savefig` on the full figure so axes, lines, legends, colorbars, annotations, and available graphics data remain editable. PNG is the preview for reports and viewing. User-specified deliverables take precedence.

PDF is optional and requires an explicit user request or an explicit requirement in the current task. Do not create PDF, SVG, EPS, or other extra formats merely because they are supported.

Keep original `.fig` files. When revising an existing figure, use a new stem or output directory unless overwriting is authorized. Do not overwrite only one member of an existing pair. Use collision-free paths for all intended outputs.

`.fig` is a MATLAB graphics-object artifact. Reopen only figures with known, trusted provenance. It is not a substitute for raw data, data transformations, or plotting source.

## Figure construction

Use explicit handles; avoid incidental `gcf`/`gca` and global `groot` styling changes. Use `tiledlayout`/`nexttile` for new multi-panel layouts where suitable. Set units, sizes, fonts, limits, colormap, colorbar labels, line styles, and markers on the owned objects. Keep comparable panels on matching normalization/color limits.

Label physical quantities and units. For this project's depth plots, follow the actual positive-down depth convention. Preserve NaN gaps, declare normalization/reference levels, and distinguish missing values from zero. Do not interpolate missing scientific results to complete a figure.

## Save the complete figure

The following assumes `fig` is an owned, completed figure, `outDir` is an explicit absolute case output directory, and `baseName` is a filename stem without separators. Resolve these in the workflow, independently of MATLAB's current directory.

```matlab
if ~isfolder(outDir)
    mkdir(outDir);
end
figPath = fullfile(outDir, [baseName '.fig']);
pngPath = fullfile(outDir, [baseName '.png']);
if isfile(figPath) || isfile(pngPath)
    error('scientific_visualization:OutputExists', ...
        'Use a new output stem or directory to preserve existing figures.');
end

drawnow;
savefig(fig, figPath);
exportgraphics(fig, pngPath, 'Resolution', 300);
```

The 300 DPI value is a general preview starting point. Choose figure dimensions and resolution for the actual destination; do not claim that a DPI value proves publication suitability.

Save after finalizing layout/annotations and before closing the figure. Preserve failures and report incomplete output pairs rather than claiming successful delivery. Close only workflow-owned handles after saving when appropriate; do not use `close all`. Do not change visibility of unrelated user figures.

In this repository, use case-specific `results/visualization/`, `results/validation/`, or another applicable results category. Keep reusable plotting code and workflow entrypoints in their established locations. Do not change the MATLAB working directory to make exports work.

## Minimal verification

For a representative newly generated figure:

1. Confirm both files exist and are nonempty, and inspect PNG dimensions and appearance.
2. Reopen the trusted `.fig` with `openfig(figPath, 'new', 'invisible')` when testing editability. Compare representative X/Y data, labels, and essential layout objects; retain NaN gaps.
3. Close only the verification copy and any owned test figure.
4. Confirm the run produced no unrequested PDF.

A plotting/save smoke test does not validate the physical model that supplied real scientific data.

## Optional PDF export

Only when PDF is a required deliverable, use the connected release's supported `exportgraphics` options and check the output for clipping and rasterized content. PDF export is not part of the default save example.
