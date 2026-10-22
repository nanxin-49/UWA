function output = sv_skill_smoke(outDir)
% Minimal artificial-data test of the skill's actual save example.
originalDirectory = pwd;
existingFigures = findall(groot, 'Type', 'figure');
x = 0:4;
y = [0 1 NaN -1 0];
fig = figure('Visible', 'off', 'Color', 'white', ...
    'Units', 'pixels', 'Position', [100 100 640 420]);
cleanupFigure = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig, 1, 1, 'Padding', 'compact');
ax = nexttile(layout);
plot(ax, x, y, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Synthetic observations', 'Tag', 'skill_smoke_line');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Amplitude (unitless)');
title(ax, 'MATLAB FIG/PNG save test');
legend(ax, 'show');
xlim(ax, [0 4]);
ylim(ax, [-1.2 1.2]);
grid(ax, 'on');
baseName = 'matlab_native';
% DOCUMENTED_SAVE_BLOCK
figInfo = dir(figPath);
pngInfo = dir(pngPath);
assert(isscalar(figInfo) && figInfo.bytes > 0);
assert(isscalar(pngInfo) && pngInfo.bytes > 0);
imageInfo = imfinfo(pngPath);
assert(imageInfo.Width > 0 && imageInfo.Height > 0);

copy = openfig(figPath, 'new', 'invisible');
cleanupCopy = onCleanup(@() close(copy)); %#ok<NASGU>
savedLine = findall(copy, 'Type', 'line', 'Tag', 'skill_smoke_line');
assert(isscalar(savedLine));
assert(isequaln(savedLine.XData, x));
assert(isequaln(savedLine.YData, y));
assert(strcmp(savedLine.Marker, 'o'));
savedAxes = ancestor(savedLine, 'axes');
assert(strcmp(savedAxes.XLabel.String, 'Time (s)'));
assert(strcmp(savedAxes.YLabel.String, 'Amplitude (unitless)'));
assert(isequal(savedAxes.XLim, [0 4]));
assert(numel(findall(copy, 'Type', 'legend')) == 1);

collisionRejected = false;
try
% DOCUMENTED_SAVE_BLOCK
catch exception
    if ~strcmp(exception.identifier, 'scientific_visualization:OutputExists')
        rethrow(exception);
    end
    collisionRejected = true;
end
assert(collisionRejected);
assert(isequal(figInfo, dir(figPath)) && isequal(pngInfo, dir(pngPath)));
assert(isempty(dir(fullfile(outDir, '*.pdf'))));
assert(strcmp(pwd, originalDirectory));
assert(all(isgraphics(existingFigures)));

output = struct('status', 'PASS', 'matlab_version', version, ...
    'fig_roundtrip_data_labels_legend', true, ...
    'nan_gap_preserved', true, 'existing_output_collision_rejected', true, ...
    'png_width', imageInfo.Width, 'png_height', imageInfo.Height, ...
    'generated_pdfs', 0, 'original_figures_preserved', true, ...
    'working_directory_preserved', true);
file = fopen(fullfile(outDir, 'matlab_smoke_result.json'), 'w');
assert(file >= 0);
cleanupFile = onCleanup(@() fclose(file)); %#ok<NASGU>
fprintf(file, '%s\n', jsonencode(output));
disp(jsonencode(output));
end
