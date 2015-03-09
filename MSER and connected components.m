colorImage = imread('137.jpg');
figure; imshow(colorImage); title('Original image')

% Detect and extract regions
grayImage = rgb2gray(colorImage);
mserRegions = detectMSERFeatures(grayImage,'RegionAreaRange',[150 2000]);
mserRegionsPixels = vertcat(cell2mat(mserRegions.PixelList));  % extract regions

% Visualize the MSER regions overlaid on the original image
figure; imshow(colorImage); hold on;
plot(mserRegions, 'showPixelList', true,'showEllipses',false);
title('MSER regions');

% Convert MSER pixel lists to a binary mask
mserMask = false(size(grayImage));
ind = sub2ind(size(mserMask), mserRegionsPixels(:,2), mserRegionsPixels(:,1));
mserMask(ind) = true;

% Run the edge detector
edgeMask = edge(grayImage, 'Canny');

% Find intersection between edges and MSER regions
edgeAndMSERIntersection = edgeMask & mserMask;
figure; imshowpair(edgeMask, edgeAndMSERIntersection, 'montage');
title('Canny edges and intersection of canny edges with MSER regions')

% [~, gDir] = imgradient(grayImage);
% You must specify if the text is light on dark background or vice versa
% gradientGrownEdgesMask = helperGrowEdges(edgeAndMSERIntersection, gDir, 'LightTextOnDark');
% figure; imshow(gradientGrownEdgesMask); title('Edges grown along gradient direction')


% Remove gradient grown edge pixels
edgeEnhancedMSERMask = ~edgeAndMSERIntersection & mserMask;

% Visualize the effect of segmentation
figure; imshowpair(mserMask, edgeEnhancedMSERMask, 'montage');
title('Original MSER regions and segmented MSER regions')

afterStrokeWidthTextMask = edgeEnhancedMSERMask;


se1=strel('disk',25);
se2=strel('disk',7);

afterMorphologyMask = imclose(afterStrokeWidthTextMask,se1);
afterMorphologyMask = imopen(afterMorphologyMask,se2);

% Display image region corresponding to afterMorphologyMask
displayImage = colorImage;
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
