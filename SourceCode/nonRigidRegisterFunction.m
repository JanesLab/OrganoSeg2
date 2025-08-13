function [outputs] = nonRigidRegisterFunction(batchStarts,batches, numWorkers, parameters, necessaryInfo)


numIterations = parameters.NumIterations;
            
smoothing = parameters.Smoothing;




regMethod = parameters.RegMethod;
bulkIFull = necessaryInfo.BulkI;
imgMaskedFull = necessaryInfo.ImgMasked;
centersFull = necessaryInfo.Centers;
ccFull = necessaryInfo.CC;
numImages = length(bulkIFull);

if length(batches) > 1
    pw = PoolWaitbar(numImages, 'Registering Images');
else
    [startInd, endInd] = getBatchEnds(batches(1), batchStarts, numImages);
    pw = PoolWaitbar(endInd-startInd+1, 'Registering Images');
end

pool = gcp("nocreate");
if isempty(pool) & numWorkers > 0
    parpool("Threads")
end
    


%h = waitbar(0,sprintf('Please wait... registering batch %d/%d', 1, length(app.BatchStarts)));
parfor (batch = batches, numWorkers)
 
    
    [startInd, endInd] = getBatchEnds(batch, batchStarts, numImages);
    
    
    
    
    
    bulkI = bulkIFull(startInd:endInd);
   
    cc = ccFull(startInd:endInd);
    
    imgMasked = imgMaskedFull(startInd:endInd);
   
    centers = centersFull(startInd:endInd);
    adjustedCenters = cell(1, endInd-startInd+1);
    displacementFields = cell(1, endInd-startInd+1);
    oldRegistered = cell(1, endInd-startInd+1);
    registeredBulkI = cell(1, endInd-startInd+1);
    registeredPerimeter = cell(1, endInd-startInd+1);


    
    %waitbar(0,h,sprintf('Please wait... registering batch %d/%d', batch, length(batchStarts)));
    adjustedCenters{1} = centers{1};
    firstImageCCLabel = labelmatrix(cc{1});
    
    %Determine whether alignment is based off of image itself or
    %based on the image mask 
    
    if strcmp(regMethod, "Mask")
        registeredBulkI{1} = logical(firstImageCCLabel);
        
        
    else
%                 background = imfilter(app.BulkI{1},fspecial('average',[500 500]),'replicate');
%                 mask = logical(firstImageCCLabel);
%                 background(mask) = 255;
%                 app.RegisteredBulkI{1} = background;
%                 newImage = int16(app.BulkI{1});
%                 mask = logical(firstImageCCLabel);
%                 %negativeImage = -newImage;
%                 %newImage(~mask) = negativeImage(~mask);
%                 newImage(mask) = 2000;
%                 app.RegisteredBulkI{1} = newImage;
        registeredBulkI{1} = bulkI{1};
        %app.RegisteredBulkI{app.BatchStarts(batch)} = imbinarize(app.BulkI{app.BatchStarts(batch)}, 'adaptive', 'ForegroundPolarity', 'Dark', 'Sensitivity', 0.6);

        
    end
    displacementFields{1} = zeros([size(bulkI{1}), 2]);
    firstPerim = bwperim(firstImageCCLabel);
    firstPerim = imdilate(firstPerim, ones(3,3));
    registeredPerimeter{1} = firstPerim;

    [meshX, meshY] = meshgrid(1:size(firstImageCCLabel,2),1:size(firstImageCCLabel,1));
    
    increment(pw)
    for selected = 2:endInd-startInd+1
        prev = selected - 1;
%                 

        %firstImageCCLabel = labelmatrix(app.CC{prev});
        secondImageCCLabel = labelmatrix(cc{selected});
        
        %Perform registration with the previously registered image
        %so that displacement fields for all images align
        
        baseImage = registeredBulkI{prev};
        
        if strcmp(regMethod, 'Mask')
            currImage = logical(secondImageCCLabel);
            
        else
%                     background = imfilter(app.BulkI{selected},fspecial('average',[500 500]),'replicate');
%                     mask = logical(firstImageCCLabel);
%                     background(mask) = 255;
%                     currImage = background;
%                     newImage = int16(app.BulkI{selected});
%                     mask = logical(secondImageCCLabel);
%                     %negativeImage = -newImage;
%                     %newImage(~mask) = negativeImage(~mask);
%                     newImage(mask) = 2000;
%                     currImage = newImage;
            currImage = bulkI{selected};
            %currImage = imbinarize(app.BulkI{selected}, 'adaptive', 'ForegroundPolarity', 'Dark', 'Sensitivity', 0.6);
        end
        
        

        
        pyramids = getMaxPyramids(currImage, baseImage);

        %Create vector for iterations, reducing iterations at
        %higher levels (which are the greatest time burden)
        %improves runtime without much effect on performance
        iterations = round(linspace(numIterations, 10, pyramids));

        %Perform non-rigid registration using the previous
        %registered image, allows transformations to combine so all
        %images align
        
        [D,currReg] = imregdemons(currImage, baseImage, iterations,...
'PyramidLevels', pyramids, 'AccumulatedFieldSmoothing', smoothing, 'DisplayWaitbar', false);

        displacementFields{selected} = D;


        registeredBulkI{selected} = currReg;
        perim = bwperim(secondImageCCLabel);                
        
        
        perimDil = imdilate(perim,ones(3,3));
        perimReg = imwarp(perimDil, D);
        %Create an overlay of the previous registered image and the
        %with the perimeter overlay to asses fit
        oldRegisteredImage = imwarp(imgMasked{prev}, displacementFields{prev});
        oldRegistered{prev} = oldRegisteredImage;
        registeredPerimeter{selected} = perimReg;
        

        %Adjust centers so all align
        centersXRounded = round(centers{selected}(:,1));
        centersYRounded = round(centers{selected}(:,2));




        adjustedX = zeros(1, length(centersXRounded));
        adjustedY = zeros(1, length(centersYRounded));
        
        sourceX = meshX(:,:) + D(:,:,1);
        sourceY = meshY(:,:) + D(:,:,2);
         
%                 for x_val = 1:size(D,2)
%                     for y_val = 1:size(D,1)
%                         sourceX(y_val, x_val) = x_val + D(y_val, x_val, 1);
%                         sourceY(y_val, x_val) = y_val + D(y_val, x_val, 2);
% 
%                     end
%                 end
        
        
        %Apply the displacement field to the centers to align
        %them for comparison
        for center = 1:length(centersXRounded)
            diffFromSource = abs(sourceX-centersXRounded(center)).^2 + abs(sourceY-centersYRounded(center)).^2;
            [~,closestIndex] = min(diffFromSource(:),[],'all',"linear");
            [y, x] = ind2sub(size(D), closestIndex);
            adjustedX(center) = x;
            adjustedY(center) = y;
%                     adjustedX(center) = centersXRounded(center) - D(centersYRounded(center), centersXRounded(center), 1);
%                     adjustedY(center) = centersYRounded(center) - D(centersYRounded(center), centersXRounded(center), 2);
            
        end


        

        
        
        %adjusted = [centersXRounded adjustedX' centersYRounded adjustedY']
        adjustedCenters{selected} = [adjustedX' adjustedY'];
        
        increment(pw)


    end
    outputs(batch).adjustedCentersFull = adjustedCenters;


    outputs(batch).registeredBulkIFull = registeredBulkI;
    outputs(batch).displacementFieldsFull = displacementFields;
    
    outputs(batch).registeredPerimeterFull = registeredPerimeter;
    
    outputs(batch).oldRegisteredFull = oldRegistered;
    
end
end


function [startInd, endInd] = getBatchEnds(batch, batchStarts, numImages)
            
    startInd = batchStarts(batch);
    
    if batch<length(batchStarts)
        endInd = batchStarts(batch+1)-1;
    else
        endInd = numImages;
    end
end


function maxPyr = getMaxPyramids(curr, base)
            dims = [size(curr), size(base)];
            minDim = min(dims);
            
            maxPyr = floor(log2(minDim));
        end