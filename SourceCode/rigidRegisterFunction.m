function [outputs] = rigidRegisterFunction(batchStarts,batches, numWorkers, parameters, necessaryInfo)




regMethod = parameters.RegMethod;
bulkIFull = necessaryInfo.BulkI;
numImages = length(bulkIFull);
imgMaskedFull = necessaryInfo.ImgMasked;
centersFull = necessaryInfo.Centers;
ccFull = necessaryInfo.CC;
registrationType = parameters.RegistrationType;
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
            


parfor (batch = batches, numWorkers)
 
    
    [startInd, endInd] = getBatchEnds(batch, batchStarts, numImages);
    startInd
    
    
    
    
    bulkI = bulkIFull(startInd:endInd);
   
    cc = ccFull(startInd:endInd);
    
    imgMasked = imgMaskedFull(startInd:endInd);
   
    centers = centersFull(startInd:endInd);
    adjustedCenters = cell(1, endInd-startInd+1);
    totalTform = cell(1, endInd-startInd+1);
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

        
        registeredBulkI{1} = bulkI{1};
       
    end
    totalTform{1} = affine2d(eye(3));
    firstPerim = bwperim(firstImageCCLabel);
    firstPerim = imdilate(firstPerim, ones(3,3));
    registeredPerimeter{1} = firstPerim;


    
    increment(pw)
    for selected = 2:endInd-startInd+1
        prev = selected - 1;
         

      
        secondImageCCLabel = labelmatrix(cc{selected});
        
        %Perform registration with the previously registered image
        %so that displacement fields for all images align
        
        baseImage = registeredBulkI{prev};
     
        if strcmp(regMethod, 'Mask')
            currImage = logical(secondImageCCLabel);
            
        else
            currImage = bulkI{selected};
        end

     
        if contains(registrationType, 'Rotate')
            [totalTformCurr, ~] = imregcorr(currImage, baseImage, 'rigid');
        else
            [totalTformCurr, ~] = imregcorr(currImage, baseImage, 'translation');
        end
      

%                 orig = tform.T
%                 totalTform = affine2d(app.TotalTform{prev}.T*tform.T);             
        totalTform{selected} = totalTformCurr;

        
        sameAsInput = affineOutputView(size(baseImage), totalTformCurr, "BoundsStyle","SameAsInput");
        registeredBulkI{selected} = imwarp(currImage, totalTform{selected},"OutputView",sameAsInput);
        perim = bwperim(secondImageCCLabel);                
        perimDil = imdilate(perim,ones(3,3));
        perimReg = imwarp(perimDil, totalTformCurr, "OutputView",sameAsInput);
     
        %Create an overlay of the previous registered image and the
        %with the perimeter overlay to asses fit
        oldRegisteredImage = imwarp(imgMasked{prev}, totalTform{prev}, "OutputView",sameAsInput);
        oldRegistered{prev} = oldRegisteredImage;
        registeredPerimeter{selected} = perimReg;
 

        %Adjust centers so all align
        centersXRounded = round(centers{selected}(:,1));
        centersYRounded = round(centers{selected}(:,2));


       

        adjustedX = zeros(1, length(centersXRounded));
        adjustedY = zeros(1, length(centersYRounded));
        
        
         
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
            oldCoord = [centersXRounded(center), centersYRounded(center), 1];
            newCoord = oldCoord*totalTformCurr.T;
            
            adjustedX(center) = newCoord(1);
            adjustedY(center) = newCoord(2);
%                     adjustedX(center) = centersXRounded(center) - D(centersYRounded(center), centersXRounded(center), 1);
%                     adjustedY(center) = centersYRounded(center) - D(centersYRounded(center), centersXRounded(center), 2);
            
        end

        %adjusted = [centersXRounded adjustedX' centersYRounded adjustedY']
        adjustedCenters{selected} = [adjustedX' adjustedY'];
        %waitbar((selected-1)/(length(bulkI)-1))
        increment(pw)
     


        

        
    end
   
    
    
    outputsCurr = struct;
    outputsCurr.adjustedCentersFull = adjustedCenters;


    outputsCurr.registeredBulkIFull = registeredBulkI;
    outputsCurr.totalTformFull = totalTform;
    
    outputsCurr.registeredPerimeterFull = registeredPerimeter;
    
    outputsCurr.oldRegisteredFull = oldRegistered;
    outputs(batch) = outputsCurr;
    
    
    
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