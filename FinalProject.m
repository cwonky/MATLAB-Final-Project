%% MATLAB FINAL PROJECT: Color by CAIT
% Names: Amanda Lo, Crystal Truong, Isabella Gaeta, Tiffany He
clear all; close all; clc;

% 1. Solicit image from user
imageName = input("Enter image file name: ",'s');
image = imread(imageName);

% Looping through channels allows imadjust to work on color photos
for c = 1:3
    image(:,:,c) = imadjust(image(:,:,c));
end

% 2. Filter image with sharp settings
smoothedImage = imbilatfilt(image, 'DegreeOfSmoothing', 0.001);
denoisedImage = medfilt3(smoothedImage, [3 3 3]);
filteredImage = imopen(denoisedImage, strel('disk', 0));

% 3. Indexing & Colormap Selection
availableColormaps = {'hot','cool','parula','turbo','jet', 'lines', 'hsv'};
choiceIdx = listdlg('ListString', availableColormaps, 'PromptString', 'Select a colormap:','SelectionMode','single');
if isempty(choiceIdx), error('No colormap selected.'); end
colormapChoice = availableColormaps{choiceIdx};

% 4. Prepare Palette
Ncolors = 20; 
try
    cmFunc = str2func(colormapChoice);
    highRes = cmFunc(1024);
catch
    highRes = parula(1024);
end
idxSample = round(linspace(1, size(highRes,1), Ncolors));
cm = highRes(idxSample, :);

% 5. Color Mapping
imgD = im2double(filteredImage);
if ndims(imgD) == 2, imgD = repmat(imgD, [1,1,3]); end
px = reshape(imgD, [], 3);
D = pdist2(px, cm, 'squaredeuclidean'); 

temp = 0.01; 
W = exp(-D / temp);
W = W ./ sum(W, 2); 
softColors = W * cm; 

[~, idxPrimary] = min(D, [], 2);
indexedMat = reshape(uint8(idxPrimary), size(imgD,1), size(imgD,2));

% 6. Post-process Regions
minRegionArea = 5;    
for k = 1:Ncolors
    mask = indexedMat == k;
    if ~any(mask(:)), continue; end
    mask = imclose(mask, strel('disk', 2));
    mask = bwareaopen(mask, minRegionArea);
    mask = imfill(mask, 'holes');
    indexedMat(mask) = k;
end
indexedMat = uint8(max(1, min(Ncolors, round(indexedMat))));

% 7. Generate Outlines
finalColoredDiscrete = im2uint8(ind2rgb(double(indexedMat), cm));
outline = false(size(indexedMat));
for k = 1:Ncolors
    mask = indexedMat == k;
    if any(mask(:)), outline = outline | bwperim(mask, 8); end
end
outline = imdilate(outline, strel('disk', 0));

% 8. Initilize canvas variables
% Create the white background before starting the number loop
whiteBg = uint8(255 * ones(size(finalColoredDiscrete)));
textImage = whiteBg; % Defines textImage so the loop can see it

% Paint outlines black on the textImage
for c = 1:3
    layer = textImage(:,:,c);
    layer(outline) = 0;
    textImage(:,:,c) = layer;
end

% 9. Insert numbers using distance transform
for k = 1:Ncolors
    mask = indexedMat == k;
    if any(mask(:))
        [L, numComp] = bwlabel(mask, 8);
        for comp = 1:numComp
            compMask = (L == comp);
            if nnz(compMask) < minRegionArea, continue; end
            
            % Find the "deepest" point for the number
            distMap = bwdist(~compMask); 
            [~, maxIdx] = max(distMap(:));
            [y, x] = ind2sub(size(compMask), maxIdx);
            
            % Render text
            pos = [x-6, y-8];
            tmpFrame = uint8(255 * ones(size(textImage)));
            tmpFrame = insertText(tmpFrame, pos, num2str(k), 'FontSize', 12, ...
                'BoxOpacity', 0, 'TextColor', 'black', 'AnchorPoint', 'LeftTop');
            
            % Apply text via grayscale mask to keep it pure black
            grayTmp = rgb2gray(tmpFrame);
            textMask = grayTmp < 250;
            for ch = 1:3
                layer = textImage(:,:,ch);
                layer(textMask) = 0;
                textImage(:,:,ch) = layer;
            end
        end
    end
end

% 10. Display results
figure('Name', 'The Canvas'); imshow(textImage);
figure('Name', 'Color Reference'); imshow(finalColoredDiscrete);

% 11. Make and display color key
hKey = figure('Name', 'Color Key'); hold on;
for k = 1:Ncolors
    rectangle('Position', [0.2, Ncolors-k, 0.6, 0.7], 'FaceColor', cm(k,:), 'EdgeColor', 'k');
    text(1.0, Ncolors-k+0.35, sprintf('Color #%d', k), 'FontSize', 12, 'FontWeight', 'bold');
end
axis([0 3 -1 Ncolors]); axis off;
