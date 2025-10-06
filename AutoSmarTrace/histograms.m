function histograms(matFile, outputDir, varargin)
% histograms(matFile, outputDir, Name,Value,...) builds angle
% histograms and Gaussian fits from SmarTrace sampled segment data.
%
% Required:
% - matFile: path to a SmarTrace re-analysis .mat (contains sampled_segments)
% - outputDir: folder to save figures and summary files
%
% Optional Name-Value pairs:
% - 'AngleRange'       [default pi]       range of angles to wrap to [-R, R]
% - 'NumBins'          [default 50]       number of histogram bins
% - 'LengthBinSize'    [default 20]       nm per length bin for per-bin plots
% - 'MinCountsPerBin'  [default 40]       minimum samples to fit a bin
%
% Examples:
%   histograms('40 Lp/SmTrAngleTest40lp301_SmTr.mat', '40 Lp/Histograms');
%   histograms('file.mat', 'out', 'NumBins', 60, 'LengthBinSize', 10);

% Tolerate a single stray positional argument (e.g., '0') by ignoring it.
if ~isempty(varargin) && mod(numel(varargin),2)==1
    warning('histograms:ignoringTrailingArg', ...
        'Ignoring trailing argument without a value: %s', mat2str(varargin{end}));
    varargin(end) = [];
end

% Dispatch to the core implementation (kept below as subfunction)
histogram_gaussian_analysis_core(matFile, outputDir, varargin{:});
end

function histogram_gaussian_analysis_core(matFile, outputDir, varargin)
% Core implementation previously named histogram_gaussian_analysis

p = inputParser;
addRequired(p, 'matFile', @(x) ischar(x) || isstring(x));
addRequired(p, 'outputDir', @(x) ischar(x) || isstring(x));
addParameter(p, 'AngleRange', pi, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NumBins', 50, @(x) isnumeric(x) && isscalar(x) && x > 4);
addParameter(p, 'LengthBinSize', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MinCountsPerBin', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parse(p, matFile, outputDir, varargin{:});

opts = p.Results;
matFile = char(matFile);
outputDir = char(outputDir);

if ~exist(matFile, 'file')
    error('histograms:FileNotFound', 'Cannot find %s', matFile);
end

data = load(matFile);
if ~isfield(data, 'sampled_segments')
    error('histograms:MissingData', 'File %s does not contain sampled_segments.', matFile);
end

segments = data.sampled_segments(:);
if isempty(segments)
    error('histograms:EmptySegments', 'sampled_segments is empty in %s.', matFile);
end

thetaVals = [];
lengthVals = [];
drawVals = [];

for idx = 1:numel(segments)
    seg = segments{idx};
    if ~isstruct(seg)
        continue;
    end
    if ~isfield(seg, 'theta') || ~isfield(seg, 'sep')
        continue;
    end
    th = double(seg.theta);
    sp = double(seg.sep);
    if ~isfinite(th) || ~isfinite(sp)
        continue;
    end
    thetaVals(end+1,1) = th; %#ok<AGROW>
    lengthVals(end+1,1) = sp; %#ok<AGROW>
    if isfield(seg, 'draw') && isfinite(double(seg.draw))
        drawVals(end+1,1) = double(seg.draw); %#ok<AGROW>
    else
        drawVals(end+1,1) = 1; %#ok<AGROW>
    end
end

if numel(thetaVals) < 5
    error('histograms:InsufficientData', 'Not enough valid segments were found in %s.', matFile);
end

thetaVals = normalizeAngles(thetaVals, opts.AngleRange);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

summary.overall = buildHistogramAndFit(thetaVals, opts.NumBins, opts.AngleRange, fullfile(outputDir, 'angle_hist_overall.png'), 'Overall');

lenEdges = 0:opts.LengthBinSize:(max(lengthVals)+opts.LengthBinSize);
lenIdx = discretize(lengthVals, lenEdges);
summary.perLength = struct('binLabel', {}, 'binEdges', {}, 'mu', {}, 'sigma', {}, 'count', {}, ...
    'histCounts', {}, 'histEdges', {}, 'binWidth', {}, 'figPath', {});

for bin = 1:max(lenIdx)
    mask = (lenIdx == bin);
    if nnz(mask) < opts.MinCountsPerBin
        continue;
    end
    binAngles = thetaVals(mask);
    if numel(binAngles) < 5
        continue;
    end
    label = sprintf('Length %d-%d nm', round(lenEdges(bin)), round(lenEdges(bin+1)));
    fileName = sprintf('angle_hist_len_%03d_%03d.png', round(lenEdges(bin)), round(lenEdges(bin+1)));
    per = buildHistogramAndFit(binAngles, opts.NumBins, opts.AngleRange, fullfile(outputDir, fileName), label);
    per.binLabel = label;
    per.binEdges = [lenEdges(bin), lenEdges(bin+1)];
    summary.perLength(end+1) = per; %#ok<AGROW>
end

summary.meta.matFile = matFile;
summary.meta.generatedOn = datestr(now, 0);
summary.meta.angleRange = opts.AngleRange;
summary.meta.numBins = opts.NumBins;
summary.meta.lengthBinSize = opts.LengthBinSize;
summary.meta.minCountsPerBin = opts.MinCountsPerBin;
summary.meta.countTotal = numel(thetaVals);

save(fullfile(outputDir, 'angle_hist_summary.mat'), 'summary');

writeSummaryTable(summary, fullfile(outputDir, 'angle_hist_summary.csv'));

fprintf('histograms: processed %d angles from %s\n', numel(thetaVals), matFile);

end

function anglesWrapped = normalizeAngles(angles, range)
if isempty(angles)
    anglesWrapped = angles;
    return;
end
anglesWrapped = mod(angles + range, 2*range) - range;
end

function result = buildHistogramAndFit(data, numBins, angleRange, outPng, label)
edges = linspace(-angleRange, angleRange, numBins+1);
centers = edges(1:end-1) + diff(edges)/2;
[counts, ~] = histcounts(data, edges);
binWidth = edges(2) - edges(1);
density = counts / (sum(counts) * binWidth);

[muHat, sigmaHat] = estimateGaussian(data);
gaussPdf = (1/(sigmaHat*sqrt(2*pi))) * exp(-0.5*((centers - muHat)/sigmaHat).^2);

f = figure('Visible', 'off');
bar(centers, density, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.65);
hold on
plot(centers, gaussPdf, 'r-', 'LineWidth', 2);
xlabel('Angle (rad)');
ylabel('Probability Density');
title(sprintf('%s (n = %d)', label, numel(data)));
legend({'Histogram', sprintf('Gaussian (mu=%.3f, sigma=%.3f)', muHat, sigmaHat)}, 'Location', 'best');
grid on

saveas(f, outPng);
close(f);

result.mu = muHat;
result.sigma = sigmaHat;
result.count = numel(data);
result.histCounts = counts;
result.histEdges = edges;
result.binWidth = binWidth;
result.figPath = outPng;
end

function [muHat, sigmaHat] = estimateGaussian(data)
data = data(:);
muHat = mean(data, 'omitnan');
sigmaHat = std(data, 'omitnan');
if ~isfinite(sigmaHat) || sigmaHat <= 0
    sigmaHat = sqrt(sum((data - muHat).^2) / max(numel(data)-1, 1));
end
if sigmaHat <= 0
    sigmaHat = eps;
end
end

function writeSummaryTable(summary, csvPath)
headers = {'Label','Count','Mu','Sigma','EdgeStart','EdgeEnd','Figure'};
rows = {};

rows(end+1,:) = {'Overall', summary.overall.count, summary.overall.mu, summary.overall.sigma, NaN, NaN, summary.overall.figPath}; %#ok<AGROW>

for k = 1:numel(summary.perLength)
    per = summary.perLength(k);
    rows(end+1,:) = {per.binLabel, per.count, per.mu, per.sigma, per.binEdges(1), per.binEdges(2), per.figPath}; %#ok<AGROW>
end

fid = fopen(csvPath, 'w');
if fid < 0
    warning('histogram_gaussian_analysis:CSVWriteFailed', 'Unable to open %s for writing.', csvPath);
    return;
end
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n', headers{:});
for r = 1:size(rows,1)
    fprintf(fid, '%s,%d,%.10g,%.10g,%.10g,%.10g,%s\n', rows{r,:});
end
fclose(fid);
end
