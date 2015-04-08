clear all; 
close all;
clc;
%a = [];
a = [2,4,5,13,18,22,23,25,29,32,33,38,41,42,44,50,51,59,65,66,67,69,75,77,78,83,86,87,90,91,92,96,97,98,99,108,109,119,120,121,123,131,136,144,145,148,150,153,155,156,157,161,163,164,169,170,175,177,178,184,193,194,198,199,203,209,213,214,215,216,219,224,225,236,239,242,243,253,259,266,276,288,294,303,307,316,322,326,341,354];
fprintf('size of a is %d\n',size(a,1));
iterator =1;
fid = fopen('C:\Users\Karthik\Documents\aspiring researchers\test_data\expt_19.txt', 'a+');
fprintf(fid, '%5s %5s %5s %5s %5s %5s \n', 'id', 'time', 'BB1', 'BB2', 'BB3', 'BB4');

for  iteration = 1:90
        %close all;
        %clc;
        fprintf('the current value of iterator is %d and a is %d\n',iterator, a(iterator));
        
        %j = a{iterator};    
        image=imread(strcat('C:\Users\Karthik\Documents\aspiring researchers\datasets\jpg\',num2str(a(iterator)),'.jpg'));
        %figure; imshow(edgeMask);title('Edge');

        %image = imsharpen(image_input);
        gray = rgb2gray(image);
        %figure; imshow(gray); title('gray');
        
        tic;

        mserRegions = detectMSERFeatures(gray,'RegionAreaRange',[150 2000]);
        mserRegionsPixels = vertcat(cell2mat(mserRegions.PixelList));  % extract regions

        % Visualize the MSER regions overlaid on the original image
        %{
        ---------figure; imshow(image); hold on;
        plot(mserRegions, 'showPixelList', true,'showEllipses',false);
        title('MSER regions');-----
        %}


        mserMask = false(size(gray));
        ind = sub2ind(size(mserMask), mserRegionsPixels(:,2), mserRegionsPixels(:,1));
        mserMask(ind) = true;

        I = histeq(gray);
        background = imopen(I,strel('disk',1));
        I2 = I - background;

        %bw of image
        bw = im2bw(I2, graythresh(I2));
        edgeMask = edge(bw, 'Canny');

        edgeAndMSERIntersection = edgeMask & mserMask;
        
        %-------figure; imshow(edgeAndMSERIntersection);
        %title('Original MSER regions and edge Enhanced MSER regions');

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

        se = strel('line',11,20);
        gradientGrownEdgesMask = imdilate(edgeAndMSERIntersection,se);
        % Remove gradient grown edge pixels
        %edgeEnhancedMSERMask = ~gradientGrownEdgesMask & I2;
        edgeEnhancedMSERMask = ~gradientGrownEdgesMask & mserMask;
        %----------------figure; imshow(edgeEnhancedMSERMask); hold on;
        connComp = bwconncomp(edgeEnhancedMSERMask); % Find connected components
        stats = regionprops(connComp,'Area','Eccentricity','Solidity');

        % Eliminate regions that do not follow common text measurements
        regionFilteredTextMask = edgeEnhancedMSERMask;

        regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Eccentricity] > .995})) = 0;
        regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Area] < 150 | [stats.Area] > 2000})) = 0;
        regionFilteredTextMask(vertcat(connComp.PixelIdxList{[stats.Solidity] < .4})) = 0;

        %------------------------figure; imshow(regionFilteredTextMask); hold on;

        %distanceImage = bwdist(~regionFilteredTextMask);
        %figure; imshow(distanceImage);


        %SWT extracted from <https://github.com/dennis1088/gaze-text-detection/blob/master/swt/src/swt.m>
        searchDirection = -1;
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
        %im = im2double(~regionFilteredTextMask);
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
        %---------------------figure, imshow(dx, []), title('Horizontal Gradient Image');
        %-----------------------figure, imshow(dy, []), title('Vertical Gradient Image');

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

        %----------------------figure, imshow(swtMap, []), title('Stroke Width Transform: First Pass');


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

        %-----------------------figure, imshow(swtMap, []), title('Stroke Width Transform: Second Pass');
        %swtMap = imwrite(swtMap, raw, jpg);
        strokeWidthImage = swtMap; %imread('raw.jpg');

        % Find remaining connected components
        connComp = bwconncomp(strokeWidthImage & regionFilteredTextMask);
        afterStrokeWidthTextMask = regionFilteredTextMask;
        for i = 1:connComp.NumObjects
            strokewidths = strokeWidthImage(connComp.PixelIdxList{i});
            % Compute normalized stroke width variation and compare to common value
            if std(strokewidths)/mean(strokewidths) > 0.35
                afterStrokeWidthTextMask(connComp.PixelIdxList{i}) = 0; % Remove from text candidates
            end
        end
        t(iterator) = toc;
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
            %rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)], 'EdgeColor','r','LineWidth',2);
            contents = [typecast(a(iterator), 'double'); t(iterator); thisBB(1); thisBB(2); thisBB(3); thisBB(4)];
            fprintf(fid, '%5d %5d %5d %5d %5d %5d\n', contents);
            %imwrite(displayImage, strcat('C:\Users\Karthik\Documents\aspiring researchers\test_data\spr_',num2str(a(iterator)),'.jpg'));
            
        end
        
        iterator = iterator +1;
        %continue
        
end
fclose(fid);
