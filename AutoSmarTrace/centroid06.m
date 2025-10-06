function [wcv, debugOut] = centroid06(pts, drw, THRESH, MAX, debugOptions)
%pattern matching algorithm
% Copyright 2016 Naghmeh Rezaei
% ******Do NOT distribute******

if nargin < 5 || isempty(debugOptions)
    debugOptions = struct();
end

DEBUG_XCORR = get_field(debugOptions, 'enabled', false);
W_PLOT_EVERY = get_field(debugOptions, 'plotEvery', 1);
SPOT_COL = get_field(debugOptions, 'spotColumn', 'argmax');
DEBUG_XCORR_MODE = lower(char(get_field(debugOptions, 'mode', 'both')));
DEBUG_XCORR_FIGBASE = get_field(debugOptions, 'figureBase', 205);
XCORR_FIG_STRIDE = get_field(debugOptions, 'figureStride', 50);
CAPTURE_DEBUG = DEBUG_XCORR && get_field(debugOptions, 'captureData', false);
CHAIN_LABEL = get_field(debugOptions, 'chainLabel', NaN);

CENTROID_FIG_BASE = get_field(debugOptions, 'centroidFigureBase', 200);
CENTROID_FIG_STRIDE = get_field(debugOptions, 'centroidFigureStride', 20);
if ~isnan(CHAIN_LABEL)
    CENTROID_FIG_BASE = CENTROID_FIG_BASE + CHAIN_LABEL * CENTROID_FIG_STRIDE;
    DEBUG_XCORR_FIGBASE = DEBUG_XCORR_FIGBASE + CHAIN_LABEL * XCORR_FIG_STRIDE;
end

debugOut = struct('settings', debugOptions);
if CAPTURE_DEBUG
    debugOut.snapshots = cell(0,1);
else
    debugOut.snapshots = {};
end

%THRESH = 80;
if THRESH == 0 && MAX == 0
    MAX = 1;
end
pts = (pts-THRESH)/(MAX-THRESH);

ll = size(pts,2);
lc = ll/2+0.5; %finds the centre point for the width
yl = size(pts,1);
yc = round(yl/2+0.5); %lc and yc centre of the search grid

%widths = 3:ll;
edgemargin = 4;
widths = 3:ll-(edgemargin+1)*2;
wcv = zeros(max(widths),ll);

ff = gcf;
if drw
    
    figure(2);
end
for w=widths
    w2=(w-1)/2; %pixels around the centre point (width)
    
    % % parabola
    % degree 4 polynomial with width w (number of pixels)
    f2 = w2^4-((1:ll)-lc).^4;
    f2i = f2<0;
    f2(f2i) = 0;
    f2 = f2/sum(f2);

%     % % parabola
%     % degree 2 polynomial with width w (number of pixels)
%     f2 = w2^2-((1:ll)-lc).^2;
%     f2i = f2<0;
%     f2(f2i) = 0;
%     f2 = f2/sum(f2);

    % negative values on the edge to penalize the background (unmatched points)
    border_width = lc-floor(w2)-edgemargin-1;
    f2i(1:border_width)=0;
    f2i(end-border_width+1:end)=0;
    %f2(f2i) = -0.5/edgemargin;
    f2(f2i) = -0.2/edgemargin;
    
    % heuristically determine the sc width, cannot pick a large width
    % because if the chain width varies it fails to recognize the chain, if
    % it is only 1 it recognizes the noise lines
    sc = 2; 
    %so3=zeros(size(pts,1),size(pts,2));
    %so3(yc-sc:yc+sc,:) = repmat(f2, sc*2+1, 1);
    so3 = repmat(f2, sc*2+1, 1) ./ (sc*2+1);

    doDebug = DEBUG_XCORR && mod(w, W_PLOT_EVERY) == 0;
    doXCorrDebug = doDebug && ~strcmpi(DEBUG_XCORR_MODE, 'centroid-only');
    doCentroidDebug = doDebug && ~strcmpi(DEBUG_XCORR_MODE, 'xcorr-only');


% 1st debugging slice======================================================
    if doCentroidDebug
        % Figure A: inputs to the correlation at this width
        figure(CENTROID_FIG_BASE + 1); clf
        if isnan(CHAIN_LABEL)
            figTitle = sprintf('centroid06: inputs @ w=%d', w);
        else
            figTitle = sprintf('centroid06: chain %d inputs @ w=%d', CHAIN_LABEL, w);
        end
        set(gcf, 'Name', figTitle);
        tiledlayout(2,3, 'TileSpacing','compact','Padding','compact');
    
        % (A1) pts (already normalized in this file)
        nexttile; imagesc(pts, [-1 1]); axis image ij; colormap gray
        title(sprintf('pts (normalized), w=%d', w));
    
        % (A2) quartic bump f2 (with edge penalties)
        nexttile; plot(f2,'LineWidth',1.2); hold on
        xline(lc,'r--','Center'); 
        yline(0,'k:');
        title('f2 (quartic bump w/ edge penalties)');
        xlabel('Perpendicular index (1..ll)'); ylabel('Amplitude');
    
        % (A3) replicated template so3 (rows = sc*2+1)
        nexttile; imagesc(so3); axis image ij
        title(sprintf('so3 = repmat(f2, %d rows)', size(so3,1)));
        xlabel('Perpendicular (cols)'); ylabel('Rows in template');
    end
% Ending 1st slice=========================================================


    
    if doXCorrDebug
        if isnan(CHAIN_LABEL)
            suffix = sprintf('w=%d', w);
        else
            suffix = sprintf('chain %d — w=%d', CHAIN_LABEL, w);
        end
        xcorrDebugCfg = struct( ...
            'enabled', true, ...
            'figureBase', DEBUG_XCORR_FIGBASE, ...
            'mode', DEBUG_XCORR_MODE, ...
            'titleSuffix', suffix);
        [v, dbgXCorr] = xcorr2normalized(pts, so3, xcorrDebugCfg);
    else
        v = xcorr2normalized(pts, so3); %image registration with a 4th polynomial filter
        dbgXCorr = [];
    end

    if CAPTURE_DEBUG
        snap = struct('width', w, 'pts', pts, 'template', so3, ...
            'correlation', v, 'aux', dbgXCorr);
        debugOut.snapshots{end+1,1} = snap; %#ok<AGROW>
    end



% Adding 2nd slice=========================================================
    if doCentroidDebug
        % Figure B: correlation map + weighting geometry
        figure(CENTROID_FIG_BASE + 2); clf
        if isnan(CHAIN_LABEL)
            figTitle = sprintf('centroid06: corr+weights @ w=%d', w);
        else
            figTitle = sprintf('centroid06: chain %d corr+weights @ w=%d', CHAIN_LABEL, w);
        end
        set(gcf, 'Name', figTitle);
        tiledlayout(2,3, 'TileSpacing','compact','Padding','compact');
    
        % (B1) normalized cross-correlation result
        nexttile; imagesc(v, [0 1]); axis image ij; colormap gray
        title('v = xcorr2normalized(pts, so3) (without centroid)');
    
        % (B2) the Gaussian weights used to collapse rows -> 1×ll
        yl = size(pts,1); yc = round(yl/2+0.5);
        gw = normpdf(1:yl, yc, 3)';          % column vector
        gw = gw / sum(gw);                   % normalize for visualization only
        nexttile; plot(gw,'LineWidth',1.2); xlim([1 yl]);
        title(sprintf('Gaussian weights along rows (yc=%d, \\sigma=3)', yc));
        xlabel('Row index (y)'); ylabel('weight');
    
        % (B3) pick one column to show the dot-product geometry
        switch SPOT_COL
            case 'argmax'
                [~, j0] = max(v(yc,:));      % column with max at center row
            case 'center'
                j0 = round(size(v,2)/2);
            otherwise
                j0 = max(1, min(size(v,2), round(SPOT_COL)));
        end
        colv = v(:, j0);
        nexttile; 
        plot(colv,'LineWidth',1.2); hold on
        yyaxis right; plot(gw,'--','LineWidth',1.2); yyaxis left
        title(sprintf('Column j=%d: v(:,j) & Gaussian weights', j0));
        xlabel('Row index (y)'); ylabel('v(:,j)');
    
        % (B4) the weighted projection across rows (this becomes wcv(w,:))
        wcv_row_dbg = (normpdf(1:yl, yc, 3) * v);   % same op as code below
        nexttile([1 2]); plot(wcv_row_dbg, 'LineWidth',1.2);
        xline(j0,'k:');
        title('Projected row = gw^T * v (with centroid)');
        xlabel('Perpendicular index (cols)'); ylabel('score');
    
        drawnow;
    end

% Ending 2nd slice=========================================================





    if size(v)~=size(pts)
        error('Size of the result does not match')
    end
    
    %almost equivalent to v(yc) but values of the neighbours are taken into
    %account. it helps to smooth out the noise
    % returns a width x ll matrix, corr values in y direction multiplied by
    %gauss sum (weighted gaussian sum)
    wcv(w,:) = normpdf(1:yl, yc, 3) * v; %at (width,centre)= score value
    
    if drw
        figure(2);
        subplot(2,3,1);
        hold off;
        imagesc(pts, [-1 1]);
        axis square;
        %hold on;
        %plot(c,yc,'r*');
        %plot(c-w2,yc,'m.');
        %plot(c+w2,yc,'m.');

        subplot(2,3,2);
        imagesc(v, [0 1]);
        %surf(v);
        axis square;
        %title(sprintf('width:%d val:%2.3f', w, m))
        title(sprintf('width:%d', w))
        
        subplot(2,3,3);
        %hold on;
        %plot3(repmat(w, 1, ll), 1:ll, v2, '.');
        imagesc(wcv, [0 1]);
        %surf(wcv);

        subplot(2,3,4);
        surf(pts);
        zlabel('Normalized Intensity')
        xlabel('Perpendicular direction (px)')
        ylabel('Tangential direction (px)')
        
%         %%to plot heat plot of a subsection of the chian:
%         figure();imagesc(pts)
%         set(gca,'dataAspectRatio',[1 1 1])
%         set(gca,'XTick',1:4:49,'XTickLabels',-6:6)
%         set(gca,'YTick',1:4:17,'YTickLabels',-2:2)
%         hold on
%         plot(lc,yc,'k*','markersize',12)
        %%%%
        
%         if w==30
%             subplot(2,3,4);
%             imagesc(v, [0 2]);
%             axis square;
%             %title(sprintf('width:%d val:%2.3f', w, m))
%             title(sprintf('width:%d', w))
%         end

%         if w==40
%             subplot(2,3,5);
%             imagesc(v, [0 1]);
%             axis square;
%             %title(sprintf('width:%d val:%2.3f', w, m))
%             title(sprintf('width:%d', w))
%         end
        
        
        subplot(2,3,6);
        surf(so3);
        set(gca,'XTick',1:4:49,'XTickLabels',-6:6)
set(gca,'YTick',1:6,'YTickLabels',0:0.25:1)

        set(gca,'fontsize',13)
xlabel('Perpendicular direction (px)')
ylabel('Tangential direction (px)')
zlabel('Polynomial Value')

%         pause(0.5);
    end
end

if drw
    figure(ff);
end
end

function value = get_field(opts, name, defaultValue)
%GET_FIELD Safe struct field accessor with default fallback.
    if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    else
        value = defaultValue;
    end
end
