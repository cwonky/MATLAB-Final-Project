%% MATLAB FINAL PROJECT
% Title: Color by CAIT
% Names: Amanda Lo, Crystal Truong, Isabella Gaeta, Tiffany He

clear all
close all
clc

% Solicit image from user
image = imread(input("Enter image file name: ",'s'));

% Filter image with three functions
smoothedImage = imbilatfilt(image);
denoisedImage = medfilt3(smoothedImage);
filteredImage = imopen(denoisedImage, strel('disk', 5));

% Indexing - Allow user to choose from a set of MATLAB colormaps
% Provide a balanced set covering warm, cool and diverse perceptual maps
availableColormaps = {'hot','cool','parula','turbo','jet'};
choiceIdx = listdlg('ListString', availableColormaps, 'PromptString', 'Select a colormap:','SelectionMode','single');
if isempty(choiceIdx)
    error('No colormap selected.');
end
colormapChoice = availableColormaps{choiceIdx};
% Ensure chosen name corresponds to an existing colormap function; if not, fall back to parula
if exist(colormapChoice, 'file') ~= 2 && ~ismember(colormapChoice, {'parula','hot','cool','turbo','jet'})
    colormapChoice = 'parula';
end

% Convert filtered image to double in [0,1]
if isinteger(filteredImage)
    imgD = im2double(filteredImage);
else
    imgD = filteredImage;
end

% Prepare a richer colormap and map each pixel with fine granularity to preserve detail
Ncolors = 20; % require at least 20 colors

% If the chosen colormap function exists, sample it densely; otherwise fall back to parula
try
    cmFunc = str2func(colormapChoice);
    % Sample a high-resolution colormap and then pick Ncolors evenly spaced entries
    highRes = cmFunc(256);
catch
    highRes = parula(256);
end
% Pick Ncolors evenly spaced colors from the high-resolution map to ensure smooth transitions
idxSample = round(linspace(1, size(highRes,1), Ncolors));
cm = highRes(idxSample, :);

% Ensure imgD is RGB in double [0,1] for distance mapping
if ndims(imgD) == 2
    imgRGB = repmat(imgD, [1,1,3]);
else
    imgRGB = imgD;
end
imgRGB = im2double(imgRGB);

% Instead of hard-quantizing each pixel immediately, perform local smoothing of the
% color assignment by using soft assignment (weighted average of nearest colors) to
% preserve structure while still limiting palette to Ncolors.
px = reshape(imgRGB, [], 3); % M x 3

% Compute squared Euclidean distance to each colormap color
D = pdist2(px, cm, 'squaredeuclidean'); % M x Ncolors

% Use softmax-like weighting on negative distances to compute a soft color projection.
% This avoids aggressive rounding while still biasing towards nearest palette colors.
% Temperature controls how soft/hard the assignments are:
temp = 0.05; % smaller -> harder assignment; adjust to preserve detail
W = exp(-D / temp);
W = W ./ sum(W, 2); % normalize per pixel

% Compute soft-colored pixels as weighted sum of colormap colors
softColors = W * cm; % M x 3

% For labeling/indexing purposes, also keep a primary index per pixel (nearest color)
[~, idxPrimary] = min(D, [], 2);
indexedMat = reshape(uint8(idxPrimary), size(imgRGB,1), size(imgRGB,2));

% Recompute finalColored from the softColors so image details remain visible
finalColored = reshape(softColors, size(imgRGB));
finalColored = im2uint8(finalColored);

% Post-process labeled regions to merge tiny components and close outlines,
% but be conservative so small important features are not lost.
minRegionArea = 100;    % reduce minimum area to preserve smaller shapes
closingRadius = 2;      % smaller closing to avoid over-merging

seClose = strel('disk', closingRadius);

for k = 1:Ncolors
    mask = indexedMat == k;
    if ~any(mask(:))
        continue
    end
    % Close very small gaps
    mask = imclose(mask, seClose);
    % Remove only tiny specks
    mask = bwareaopen(mask, minRegionArea);
    % Fill holes to make regions coherent
    mask = imfill(mask, 'holes');
    indexedMat(mask) = k;
end

% Ensure indices are valid and contiguous
indexedMat = uint8(max(1, min(Ncolors, round(indexedMat))));

% Recreate a discrete-colored image (for outlines/labels) using the chosen palette
% Use ind2rgb expecting indices 1..Ncolors mapped to cm
finalColoredDiscrete = ind2rgb(double(indexedMat), cm);
finalColoredDiscrete = im2uint8(finalColoredDiscrete);

% Use the soft-colored image as the main color output but keep a discrete copy for labeling/outline reference
finalColored = im2uint8(finalColored);
finalColored = imresize(finalColored, [size(imgRGB,1), size(imgRGB,2)]); %#ok<NASGU>
finalColored = finalColored; % ensure variable exists for downstream code

% Also make sure a version named finalColored (used later) exists and matches discrete colors for outlines
% Choose to prefer the soft image for display but keep discrete for indexed operations
finalColored = finalColored; %#ok<NASGU>
finalColored = finalColoredDiscrete; % use discrete-colored for outline consistency

% Take original image and make all white except for black outlines of color regions
% Create outlines from the indexed color regions, but ensure boundaries align with
% the original image intensity edges to better follow visual structure.
outline = false(size(indexedMat));

% Compute per-color outlines and accumulate
for k = 1:Ncolors
    mask = indexedMat == k;
    if any(mask(:))
        b = bwperim(mask, 8);
        outline = outline | b;
    end
end

% Close tiny gaps and ensure outlines are continuous
outline = imclose(outline, strel('disk', 1));

% Refine outlines by preferring strong gradients from the original image where available
if ndims(image) == 3
    grayOrig = rgb2gray(image);
else
    grayOrig = image;
end
grayOrig = im2double(grayOrig);

% Detect strong edges; use moderate smoothing by specifying sigma via imgaussfilt if noisy
edges = edge(grayOrig, 'Canny');

% Keep outline pixels that align with image edges, but don't discard outline where edges absent
outlineDilated = imdilate(outline, strel('disk', 1));
refined = (outline & edges) | (outline & ~edges);

% If refinement removed too many outline pixels, revert to original outline
if nnz(refined) < max(10, round(0.25 * nnz(outline)))
    refined = outline;
end
outline = bwmorph(refined, 'thin', Inf);
outline = bwareaopen(outline, 10);
outline = logical(outline);

% Build white background and paint outlines black
whiteBg = uint8(255 * ones(size(finalColored)));
for c = 1:3
    ch = whiteBg(:,:,c);
    ch(outline) = 0;
    whiteBg(:,:,c) = ch;
end

% Insert numeric labels onto the white background at centroids of sufficiently large components
whiteBgWithNumbers = whiteBg;
for k = 1:Ncolors
    mask = indexedMat == k;
    if ~any(mask(:))
        continue
    end
    [L, numComp] = bwlabel(mask, 8);
    for comp = 1:numComp
        compMask = (L == comp);
        if nnz(compMask) < minRegionArea
            continue
        end
        s = regionprops(compMask, 'Centroid');
        cxy = round(s.Centroid);
        pos = [cxy(1)-8, cxy(2)-10];
        pos(1) = max(1, min(size(whiteBgWithNumbers,2)-16, pos(1)));
        pos(2) = max(1, min(size(whiteBgWithNumbers,1)-20, pos(2)));
        txt = num2str(k);
        textColor = round(255 * cm(k, :));
        whiteBgWithNumbers = insertText(whiteBgWithNumbers, pos, txt, ...
            'FontSize', 18, 'BoxColor', 'white', 'TextColor', textColor, 'BoxOpacity', 0.6);
    end
end

whiteBg = im2uint8(whiteBgWithNumbers);

% Create the colored image that also retains the blpeack outlines
% Start from finalColored and paint outlines black
finalColoredWithOutlines = finalColored;
for c = 1:3
    ch = finalColoredWithOutlines(:,:,c);
    ch(outline) = 0;
    finalColoredWithOutlines(:,:,c) = ch;
end
finalColoredWithOutlines = im2uint8(finalColoredWithOutlines);

% Create labeled overlay: place numbers at centroids of each color region
textImage = whiteBg;
labelFontSize = 14; % reduced from 18 to make numbers smaller
boxOffset = [ -6, -8 ]; % adjust offsets for smaller font (x, y)
textBoxMargin = [12, 16]; % approximate box width/height to keep text inside image

% One label per distinct color region (connected component of same index)
% Use a pure black text drawn directly to the white image to avoid colored shadows
for k = 1:Ncolors
    mask = indexedMat == k;
    if any(mask(:))
        [L, numComp] = bwlabel(mask, 8);
        for comp = 1:numComp
            compMask = (L == comp);
            % ignore very small components (they were mostly merged already)
            if nnz(compMask) < minRegionArea
                continue
            end
            s = regionprops(compMask, 'Centroid');
            cxy = round(s.Centroid);
            % position text centered on centroid, adjusted for smaller font
            pos = [cxy(1) + boxOffset(1), cxy(2) + boxOffset(2)];
            pos(1) = max(1, min(size(textImage,2)-textBoxMargin(1), pos(1)));
            pos(2) = max(1, min(size(textImage,1)-textBoxMargin(2), pos(2)));
            txt = num2str(k);

            % Instead of insertText (which can introduce colored anti-aliased shadows),
            % render text into a temporary grayscale mask and composite it as solid black.
            % Create an RGB frame for text rendering on white background
            tmpFrame = uint8(255 * ones(size(textImage)));
            tmpFrame = insertText(tmpFrame, pos, txt, 'FontSize', labelFontSize, ...
                'BoxOpacity', 0, 'TextColor', 'black', 'AnchorPoint', 'LeftTop');

            % Convert both images to grayscale masks to find where text was drawn
            grayTmp = rgb2gray(tmpFrame);
            textMask = grayTmp < 250; % threshold to capture black text while avoiding faint colored edges

            % Paint pure black where the mask indicates text
            for ch = 1:3
                layer = textImage(:,:,ch);
                layer(textMask) = 0;
                textImage(:,:,ch) = layer;
            end
        end
    end
end
textImage = im2uint8(textImage);

% Also produce a labeled colored image (colored background with black outlines and white numbers)
labeledColored = finalColoredWithOutlines;
for k = 1:Ncolors
    mask = indexedMat == k;
    if any(mask(:))
        [L, numComp] = bwlabel(mask, 8);
        for comp = 1:numComp
            compMask = (L == comp);
            if nnz(compMask) < minRegionArea
                continue
            end
            s = regionprops(compMask, 'Centroid');
            cxy = round(s.Centroid);
            pos = [cxy(1) + boxOffset(1), cxy(2) + boxOffset(2)];
            pos(1) = max(1, min(size(labeledColored,2)-textBoxMargin(1), pos(1)));
            pos(2) = max(1, min(size(labeledColored,1)-textBoxMargin(2), pos(2)));
            txt = num2str(k);
            % For the colored image, use white text for visibility
            labeledColored = insertText(labeledColored, pos, txt, 'FontSize', labelFontSize, 'BoxOpacity', 0, 'TextColor', 'white');
        end
    end
end
labeledColored = im2uint8(labeledColored);

% Display the required final products:
% 1) whitened image with black outlines and numbers within each outlined region
figure;
imshow(textImage);
title('White Background with Black Outlines and Region Numbers');

% 2) final colored image with black outlines and white numbers
figure;
imshow(labeledColored);
title('Final Colored Image with Black Outlines and Region Numbers');

% Also display plain white outlines-only image and colored-without-numbers if desired
% (kept minimal per requirements)
