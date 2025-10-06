function [c, debug] = xcorr2normalized(x, y, debugCfg)
%XCORR2NORMALIZED Normalized 2-D cross-correlation without toolbox reliance.
%
%   C = XCORR2NORMALIZED(X, Y) computes the normalized correlation map of
%   template Y against image patch X, matching the output shape (size(X)) of
%   MATLAB's NORMXCORR2 cropped to the valid region. The implementation
%   avoids the Image Processing Toolbox so it can be used in standalone
%   deployments.
%
%   [C, DEBUG] = XCORR2NORMALIZED(X, Y) also returns intermediate products
%   useful for debugging (raw correlation, local statistics, etc.).
%
%   ... = XCORR2NORMALIZED(X, Y, DEBUGCFG) enables optional visualization.
%   Pass a struct with field DEBUGCFG.enabled=true along with optional
%   fields:
%       figureBase  - figure id for plots (default 250)
%       mode        - 'raw', 'normalized', or 'both' (default 'both')
%       titleSuffix - string appended to figure name
%
%   This mirrors the behaviour we previously relied on when prototyping
%   centroid visualisations.

    if nargin < 3 || isempty(debugCfg)
        debugCfg = struct();
    end
    if ~isfield(debugCfg, 'mode') || isempty(debugCfg.mode)
        debugCfg.mode = 'both';
    end

    wantDebugOut = nargout > 1;
    debug = struct();

    % Ensure double precision for the convolutions
    x = double(x);
    y = double(y);

    templateMean = mean(y(:));
    templateZeroMean = y - templateMean;
    templateEnergy = sum(templateZeroMean(:) .^ 2);
    templateEnergy = max(templateEnergy, eps);
    templateStd = sqrt(templateEnergy);

    % Raw correlation used for visual comparisons
    rawCorr = conv2(x, rot90(y, 2), 'same');

    % Numerator uses zero-mean template (window mean handled separately)
    numerator = conv2(x, rot90(templateZeroMean, 2), 'same');

    % Local statistics for each window
    kernel = ones(size(y));
    sum_x = conv2(x, kernel, 'same');
    sum_x2 = conv2(x .^ 2, kernel, 'same');

    n = numel(y);
    windowVar = sum_x2 - (sum_x .^ 2) / n;
    windowVar = max(windowVar, 0);
    windowStd = sqrt(windowVar);

    denom = windowStd * templateStd;
    zeroMask = denom < eps;
    denom(zeroMask) = 1; % avoid division-by-zero warnings
    c = numerator ./ denom;
    c(zeroMask) = 0;

    % Clamp to numerical bounds of normalised correlation
    c = max(min(c, 1), -1);

    if wantDebugOut || (isfield(debugCfg, 'enabled') && debugCfg.enabled)
        debug.rawCorr = rawCorr;
        debug.numerator = numerator;
        debug.denominator = denom;
        debug.windowStd = windowStd;
        debug.templateStd = templateStd;
        debug.templateZeroMean = templateZeroMean;
        debug.result = c;
        debug.kernelSize = size(y);
        debug.templateMean = templateMean;

        % Compare to MATLAB's normxcorr2 when available for reference
        if exist('normxcorr2', 'file') == 2
            builtinFull = normxcorr2(y, x);
            d1 = floor((size(builtinFull) - size(x)) / 2);
            d2 = ceil((size(builtinFull) - size(x)) / 2);
            builtinCropped = builtinFull( ...
                d1(1) + 1:end - d2(1), ...
                d1(2) + 1:end - d2(2));
            debug.normxcorr2 = builtinCropped;
        else
            debug.normxcorr2 = [];
        end
    end

    if isfield(debugCfg, 'enabled') && debugCfg.enabled
        figId = 250;
        if isfield(debugCfg, 'figureBase') && ~isempty(debugCfg.figureBase)
            figId = debugCfg.figureBase;
        end
        titleSuffix = '';
        if isfield(debugCfg, 'titleSuffix') && ~isempty(debugCfg.titleSuffix)
            titleSuffix = [' — ' char(debugCfg.titleSuffix)];
        end

        figure(figId); clf
        set(gcf, 'Name', sprintf('xcorr2normalized%s', titleSuffix));
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        modeValue = lower(char(debugCfg.mode));
        showRaw = strcmp(modeValue, 'raw') || strcmp(modeValue, 'both');
        showNorm = strcmp(modeValue, 'normalized') || strcmp(modeValue, 'both');

        nexttile
        imagesc(x);
        axis image ij
        colormap(gca, 'gray');
        title('Input patch X');
        colorbar

        nexttile
        imagesc(y);
        axis image ij
        colormap(gca, 'gray');
        title('Template Y');
        colorbar

        if showRaw
            nexttile
            imagesc(rawCorr);
            axis image ij
            colormap(gca, 'parula');
            title('Raw cross-corr (un-normalised)');
            colorbar
        else
            nexttile
            axis off
        end

        if showNorm
            nexttile
            imagesc(c, [-1, 1]);
            axis image ij
            colormap(gca, 'parula');
            title('Normalised correlation');
            colorbar
        else
            nexttile
            axis off
        end

        if isfield(debug, 'normxcorr2') && ~isempty(debug.normxcorr2)
            figure(figId + 1); clf
            set(gcf, 'Name', sprintf('xcorr2normalized vs normxcorr2%s', titleSuffix));
            tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            nexttile
            imagesc(c, [-1, 1]);
            axis image ij
            colormap(gca, 'parula');
            title('Custom implementation');
            colorbar

            nexttile
            imagesc(debug.normxcorr2, [-1, 1]);
            axis image ij
            colormap(gca, 'parula');
            title('MATLAB normxcorr2');
            colorbar

            nexttile
            diffMap = debug.normxcorr2 - c;
            maxDiff = max(abs(diffMap(:)));
            maxDiff = max(maxDiff, eps);
            imagesc(diffMap, [-maxDiff, maxDiff]);
            axis image ij
            colormap(gca, 'parula');
            title('Difference (builtin - custom)');
            colorbar
        end
        drawnow;
    end
end
