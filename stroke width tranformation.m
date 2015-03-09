clear all; 
close all;
clc;
image=imread('137.jpg');
%figure; imshow(edgeMask);title('Edge');

gray = rgb2gray(image);
%figure; imshow(gray); title('gray');
mserRegions = detectMSERFeatures(gray,'RegionAreaRange',[150 2000]);
mserRegionsPixels = vertcat(cell2mat(mserRegions.PixelList));  % extract regions
mserMask = false(size(gray));
ind = sub2ind(size(mserMask), mserRegionsPixels(:,2), mserRegionsPixels(:,1));
mserMask(ind) = true;
edgeMask = edge(gray, 'Canny');
edgeEnhancedMSERMask = edgeMask & mserMask;
figure; imshowpair(mserMask, edgeEnhancedMSERMask, 'montage');
title('Original MSER regions and edge Enhanced MSER regions');

%{
connComp = bwconncomp(edgeEnhancedMSERMask); % Find connected components
stats = regionprops(connComp,'Area','Eccentricity','Solidity');

% Eliminate regions that do not follow common text measurements
regionFilteredTextMask = edgeEnhancedMSERMask;

regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Eccentricity] > .995})) = 0;
regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Area] < 150 | [stats.Area] > 2000})) = 0;
regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Solidity] < .4})) = 0;

% Visualize results of filtering
figure; imshowpair(edgeEnhancedMSERMask, regionFilteredTextMask, 'montage');
title('Text candidates before and after region filtering')
%}

distanceImage = bwdist(~edgeEnhancedMSERMask);
figure; imshow(distanceImage);
%SWT extracted from <https://github.com/dennis1088/gaze-text-detection/blob/master/swt/src/swt.m>
searchDirection = 1;
%function [ swtMap ] = swt( image, searchDirection )
%swt Preforms stoke width transform on input image
%   A novel image operator that seeks to find the value of stroke width
%   for each image pixel.  It's use is meant for the task of text
%   detection in natural images.
%
%   im = RGB input image of size m x n x 3
%   searchDirection = gradient direction is either 1 to detect dark text on light
%   background or -1 to detect light text on dark background.
%
%   swtMap = resulting mapping of stroke withs for image pixels

% Convert image to gray scale
im = im2double(rgb2gray(image));
%im = typecast(distanceImage, 'double');
%figure, imshow(im), title('Black and White Image');

% Find edges using canny edge dector
edgeMap = edge(im, 'canny');
%figure, imshow(edgeMap), title('Edges Using Canny');

% Get all edge pixel postitions
[edgePointRows, edgePointCols] = find(edgeMap);

% Find gradient horizontal and vertical gradient
sobelMask = fspecial('sobel');
dx = imfilter(im,sobelMask);
dy = imfilter(im,sobelMask');
figure, imshow(dx, []), title('Horizontal Gradient Image');
figure, imshow(dy, []), title('Vertical Gradient Image');

% Initializing matrix of gradient direction
theta = zeros(size(edgeMap,1),size(edgeMap,2));

% Calculating theta, gradient direction, for each pixel on the image.
% ***This can be optimized by using edgePointCols and edgePointRows
% instead.***
for i=1:size(edgeMap,1)
    for j=1:size(edgeMap,2)
        if edgeMap(i,j) == 1
            tan_y = dy(i,j);
            tan_x = dx(i,j);
            theta(i,j) = atan2(tan_y, tan_x);
        end
    end
end

% Getting size of the image
[m,n] = size(edgeMap);

% Initializing Stoke Width array with infinity
swtMap = zeros(m,n);
for i=1:m
    for j=1:n
        swtMap(i,j) = inf;
    end
end

% Set the maximum stroke width, this number is variable for now but must be
% made to be more dynamic in the future
maxStrokeWidth = 350;

% Initialize container for all stoke points found
strokePointsX = zeros(size(edgePointCols));
strokePointsY = zeros(size(strokePointsX));
sizeOfStrokePoints = 0;

% Iterate through all edge points and compute stoke widths
for i=1:size(edgePointRows)
    step = 1;
    initialX = edgePointRows(i);
    initialY = edgePointCols(i);
    isStroke = 0;
    initialTheta = theta(initialX,initialY);
    sizeOfRay = 0;
    pointOfRayX = zeros(maxStrokeWidth,1);
    pointOfRayY = zeros(maxStrokeWidth,1);
    
    % Record first point of the ray
    pointOfRayX(sizeOfRay+1) = initialX;
    pointOfRayY(sizeOfRay+1) = initialY;
    
    % Increase the size of the ray
    sizeOfRay = sizeOfRay + 1;
    
    % Follow the ray
    while step < maxStrokeWidth
        nextX = round(initialX + cos(initialTheta) * searchDirection * step);
        nextY = round(initialY + sin(initialTheta) * searchDirection * step);
        
        step = step + 1;
        
        % Break loop if out of bounds.  For some reason this is really
        % slow.
        if nextX < 1 || nextY < 1 || nextX > m || nextY > n
            break
        end
        
        % Record next point of the ray
        pointOfRayX(sizeOfRay+1) = nextX;
        pointOfRayY(sizeOfRay+1) = nextY;
        
        % Increase size of the ray
        sizeOfRay = sizeOfRay + 1;
        
        % Another edge pixel has been found
        if edgeMap(nextX,nextY)
            
            oppositeTheta = theta(nextX,nextY);
            
            % Gradient direction roughtly opposite
            if abs(abs(initialTheta - oppositeTheta) - pi) < pi/2
                isStroke = 1;
                strokePointsX(sizeOfStrokePoints+1) = initialX;
                strokePointsY(sizeOfStrokePoints+1) = initialY;
                sizeOfStrokePoints = sizeOfStrokePoints + 1;
            end
            
            break
        end
    end
    
    % Edge pixel is part of stroke
    if isStroke
        
        % Calculate stoke width
        strokeWidth = sqrt((nextX - initialX)^2 + (nextY - initialY)^2);
        
        % Iterate all ray points and populate with the minimum stroke width
        for j=1:sizeOfRay
            swtMap(pointOfRayX(j),pointOfRayY(j)) = min(swtMap(pointOfRayX(j),pointOfRayY(j)),strokeWidth);
        end
    end
end

figure, imshow(swtMap, []), title('Stroke Width Transform: First Pass');

afterStrokeWidthTextMask_1 = swtMap;

%swtMap = distanceImage;

% Iterate through all stoke points for a refinement pass.  Refer to figure
% 4b in the paper.
for i=1:sizeOfStrokePoints
    step = 1;
    initialX = strokePointsX(i);
    initialY = strokePointsY(i);
    initialTheta = theta(initialX,initialY);
    sizeOfRay = 0;
    pointOfRayX = zeros(maxStrokeWidth,1);
    pointOfRayY = zeros(maxStrokeWidth,1);
    swtValues = zeros(maxStrokeWidth,1);
    sizeOfSWTValues = 0;
    
    % Record first point of the ray
    pointOfRayX(sizeOfRay+1) = initialX;
    pointOfRayY(sizeOfRay+1) = initialY;
    
    % Increase the size of the ray
    sizeOfRay = sizeOfRay + 1;
    
    % Record the swt value of first stoke point
    swtValues(sizeOfSWTValues+1) = swtMap(initialX,initialY);
    sizeOfSWTValues = sizeOfSWTValues + 1;
    
    % Follow the ray
    while step < maxStrokeWidth
        nextX = round(initialX + cos(initialTheta) * searchDirection * step);
        nextY = round(initialY + sin(initialTheta) * searchDirection * step);
        
        step = step + 1;
        
        % Record next point of the ray
        pointOfRayX(sizeOfRay+1) = nextX;
        pointOfRayY(sizeOfRay+1) = nextY;
        
        % Increase size of the ray
        sizeOfRay = sizeOfRay + 1;
        
        % Record the swt value of next stoke point
        swtValues(sizeOfSWTValues+1) = swtMap(nextX,nextY);
        sizeOfSWTValues = sizeOfSWTValues + 1;
        
        % Another edge pixel has been found
        if edgeMap(nextX,nextY)
            break
        end
    end
    
    % Calculate stoke width as the median value of all swtValues seen.
    strokeWidth = median(swtValues(1:sizeOfSWTValues));
    
    % Iterate all ray points and populate with the minimum stroke width
    for j=1:sizeOfRay
        swtMap(pointOfRayX(j),pointOfRayY(j)) = min(swtMap(pointOfRayX(j),pointOfRayY(j)),strokeWidth);
    end
    
end

figure, imshow(swtMap, []), title('Stroke Width Transform: Second Pass');
%swtMap = imwrite(swtMap, raw, jpg);
afterStrokeWidthTextMask = swtMap; %imread('raw.jpg');
%figure, imshowpair(swtMap, afterStrokeWidthTextMask, 'montage'), title('stroke width writing + distanceImage and stroke width writing')
edgeMask_1 = edge(swtMap, 'Canny');
%end
%figure, imshow(edgeMask_1, []), title('Stroke Width Transform: Second Pass and edge growing');

connComp = bwconncomp(afterStrokeWidthTextMask_1); % Find connected components
stats = regionprops(connComp,'Area','Eccentricity','Solidity');

% Eliminate regions that do not follow common text measurements
regionFilteredTextMask = edgeEnhancedMSERMask;

regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Eccentricity] > .995})) = 0;
regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Area] < 150 | [stats.Area] > 2000})) = 0;
regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Solidity] < .4})) = 0;
%----------------------

strokeWidthImage = afterStrokeWidthTextMask; % Compute stroke width image

% Show stroke width image
figure; imshow(strokeWidthImage);




% Find remaining connected components
connComp = bwconncomp(regionFilteredTextMask);
afterStrokeWidthTextMask = regionFilteredTextMask;
for i = 1:connComp.NumObjects
    strokewidths = strokeWidthImage(connComp.PixelIdxList{i});
    % Compute normalized stroke width variation and compare to common value
    if std(strokewidths)/mean(strokewidths) > 0.35
        afterStrokeWidthTextMask(connComp.PixelIdxList{i}) = 0; % Remove from text candidates
    end
end

% Visualize the effect of stroke width filtering
figure; imshowpair(regionFilteredTextMask, afterStrokeWidthTextMask,'montage');
title('Text candidates before and after stroke width filtering')
%----------------------
%disp(min(swtMap(:)));
%disp(max(swtMap(:)));
%level = graythresh(swtMap);
%disp(level);
%afterStrokeWidthTextMask = im2bw(swtMap, 1);
figure, imshow(afterStrokeWidthTextMask, []), title('swtMap to bw');
% Visualize the effect of stroke width filtering
%figure; imshowpair(swtMap, afterStrokeWidthTextMask, 'montage');
%title('stroke width writing and edge+swt')


%se1= strel('line',10,10);
se1= strel('disk',25);
se2= strel('disk',7);
%se2= strel('line',5,2);

afterMorphologyMask = imclose(afterStrokeWidthTextMask,se1);
afterMorphologyMask = imopen(afterMorphologyMask,se2);

% Display image region corresponding to afterMorphologyMask
displayImage = image;
displayImage(~repmat(afterMorphologyMask,1,1)) = 0;
figure; imshow(displayImage); title('Image region under mask created by joining individual characters')


areaThreshold = 5000; % threshold in pixels
connComp = bwconncomp(afterMorphologyMask);
stats = regionprops(connComp,'BoundingBox','Area');
boxes = round(vertcat(stats(vertcat(stats.Area) > areaThreshold).BoundingBox));
for k = 1 : length(stats)
    thisBB = stats(k).BoundingBox; 
    rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)], 'EdgeColor','r','LineWidth',2); 
end
%{
for i=1:size(boxes,1)
    figure;
    imshow(imcrop(image, boxes(i,:))); % Display segmented text
    title('Text region')
end
%}





%{
thresholdValue = graythresh(I);  %# Compute an appropriate threshold
%print('-depsc',level);
BW = gray > thresholdValue; 
imshow(BW);
%BW=binary_image(image); 
BW = ~BW; 
st = regionprops(BW, 'BoundingBox' ); 
for k = 1 : length(st)
    thisBB = st(k).BoundingBox; 
    rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)], 'EdgeColor','r','LineWidth',2); 
end
%}

%}
