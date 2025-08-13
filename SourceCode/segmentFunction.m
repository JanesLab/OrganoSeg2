function [outputs] = segmentFunction(finish, numWorkers, mainParameters, additionalParameters, bulkI)
%h = waitbar(0,'Please wait... segmenting...');
            
            %timesSegmentNewMinParam = zeros(1, length(bulkIRef));
            edgePreference = additionalParameters.edgePreference;
            useReconstructedForGradient = additionalParameters.useReconstructedForGradient;
            useInverted = additionalParameters.useInverted;
            clearBorderSplit = additionalParameters.clearBorderSplit;
            windowSizeJump = additionalParameters.windowSizeJump;
            imageReconstructionRadius = additionalParameters.imageReconstructionRadius;
            imageClosingDegree = additionalParameters.imageClosingDegree;
            watershedFilterSigma = additionalParameters.watershedFilterSigma;


            ws = mainParameters.ws;
            intensityThresh = mainParameters.intensityThresh;
            sizeThresh = mainParameters.sizeThresh;
            checkOOF = mainParameters.checkOOF;
            checkDIC = mainParameters.checkDIC;
            checkWater = mainParameters.checkWater;
            checkEdge = mainParameters.checkEdge;
            edgeIntensity = mainParameters.edgeIntensity;          
            minimumCircularity = mainParameters.minimumCircularity;
            contaminantIntensity = mainParameters.contaminantIntensity;

            
            
            


            

             

                          
            pw = PoolWaitbar(finish, 'Segmenting Images');


            
            
            
            pool = gcp("nocreate");
            if isempty(pool) & numWorkers > 0
                parpool("Threads")
            end
            parfor (index = 1:finish, numWorkers)
                
          
                
         
                
    

                
                 

               
               
                
             
               
                windowSize = [round(ws/10)*10, round(ws)/10*10];
               
                %get segmenting image%%
                origI = bulkI{index};

                %Create structure element
                se = strel('disk', imageReconstructionRadius);

                Iobrcbr = prepareImage(origI, se, useInverted);

                %Threshold
                bwEdge = zeros(size(origI,1),size(origI,2));
                fudgeFactor = intensityThresh;

                %Multi-scale checkbox?
                windowSizes = 20:windowSizeJump:windowSize(1);
                if windowSizes(end) ~= windowSize(1)
                    windowSizes = [windowSizes windowSize(1)];
                end
                
                if checkOOF == 1
                    for step = 1:length(windowSizes)
                        %original segmentation 
                        j = windowSizes(step);
                        bwEdgeIter = adaptivethreshold(Iobrcbr,[j j],fudgeFactor);
                        bwEdge = bwEdge+bwEdgeIter;
                    end
                else
                    bwEdge = adaptivethreshold(Iobrcbr,windowSize,fudgeFactor);
                    %DIC Checkbox?
                    if checkDIC == 1
                        
                        %Segment on inverted image to capture shadow vs light
                        bwEdge2 = adaptivethreshold(Iobrcbr,windowSize,fudgeFactor);

                        %Sum together
                        bwEdge = bwEdge|bwEdge2;
                        bwEdge = imclose(bwEdge,strel('disk',5));
                    end
                end

                %Remove small noise
                bwsmall = bwareaopen(bwEdge, sizeThresh);
                %Smooth
                bwsmooth = imclose(bwsmall,ones(imageClosingDegree));
                %fill holes
                bw_holes = bwfill(bwsmooth,'holes');
                bw_holes = logical(bw_holes);
                

                 %Apply edge correction by using gradient to find edges
                if checkEdge == 1
                    if strcmp(edgePreference, 'Gradient Only')

                        bw_holes = gradientOnlyEdgeCorrect(bw_holes, edgeIntensity, Iobrcbr, origI, imageReconstructionRadius, useReconstructedForGradient, useInverted,sizeThresh);
                        
                    elseif strcmp(edgePreference, 'Gradient Preserve Boundary')

                        
                        bw_holes = gradientPreserveBoundaryEdgeCorrect(bw_holes, edgeIntensity, Iobrcbr, origI, imageReconstructionRadius, useReconstructedForGradient, useInverted,sizeThresh);
                        
                        
                    %Use watershed thresholding based on original
                    %segmentation
                    elseif strcmp(edgePreference, 'Gradient + Watershed')

                        bw_holes = watershedEdgeCorrect(bw_holes, origI, imageReconstructionRadius, useReconstructedForGradient, useInverted,sizeThresh);

                    end
                end

                bw_holes = logical(bw_holes);
                
               
               
                
               
    
                %Watershed checkbox?
                
                if checkWater == 1
                    if clearBorderSplit
                        bw_holes = imclearborder(bw_holes);
                    end
                   
                    bwfinal = watershedSplit(bw_holes, watershedFilterSigma, 1, sizeThresh);                    
                else
                    %Clear border
                    bw_nobord = imclearborder(bw_holes);
                    %make logical
                    bwfinal = logical(bw_nobord);                    
                    
                    
                    
                end
                
               
 

                %Add color labels
                cc = bwconncomp(bwfinal);
                   
              
            

                necessaryData = struct;
                necessaryData.cc = cc;
                necessaryData.bulkI = origI;
                
                necessaryData.minimumCircularity = minimumCircularity;
                necessaryData.contaminantIntensity = contaminantIntensity;
                
                
                
                necessaryData.checkWater = checkWater;
                
                necessaryData.bw_nowater = bw_holes;
                
               
                outputsCurr = redisplayWithChangeToSpheroidData(necessaryData);
                outputsCurr.cc = cc;
                outputsCurr.bw_noWater = bw_holes;
                
                outputsCurr.bwfinal = bwfinal;
                outputs(index) = outputsCurr;

                % imgMaskedFull{index} = outputs.ImgMasked;
                % nonCircCCFull{index} = outputs.NonCircularCC;
                % contaminantCCFull{index} = outputs.ContaminantCC;
                % nonContaminantBooleanFull{index} = outputs.NonContaminantBoolean;
                % circularBooleanFull{index} = outputs.CircularBoolean;
                % coloredLabelsFull{index} = outputs.ColoredLabels;
                % centersFull{index} = outputs.Centers;
                % boundingBoxFull{index} = outputs.BoundingBox;
                % spheroids_orderFull{index} = outputs.Spheroids_order;
                % bwfinalFull{index} = outputs.Bwfinal;

                

                
                %timesSegmentNewMinParam(index) = toc(timerVal);
                
                increment(pw)
                
                
            end
            
          

            
end



function Iobrcbr = prepareImage(origI, se, useInverted)
            if useInverted
                origI = imcomplement(origI);
            end
            
            
            
            %Sharpen
            K = imsharpen(origI);
            %Smooth background
            K = medfilt2(K);

            %Opening and closing by reconstruction

            

            %Erode and open
            Ie = imerode(K, se);
            Iobr = imreconstruct(Ie, K);
            
            %Opening-closing by reconstruction(Iobrcbr)
            Iobrd = imdilate(Iobr, se);
            Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
            Iobrcbr = imcomplement(Iobrcbr);

            %Remove border
            Iobrcbr(1:end,1) = Iobrcbr([2 2:end],2);
            Iobrcbr(1:end,end) = Iobrcbr([2 2:end],end-1);
            Iobrcbr(1,1:end) = Iobrcbr(2,[2 2:end]);
            Iobrcbr(end,1:end) = Iobrcbr(end-1,[2 2:end]);
            
end


function outputs = redisplayWithChangeToSpheroidData(necessaryData)
            cc = necessaryData.cc;
            CClabel = labelmatrix(cc);
            %Convert CClabel to bw                       
            bwfinal = logical(CClabel);            
            


            CC_temp = findCCTemp(cc);
            centers = findCentersFromCC(cc);            
            bbox = findBoundingBoxFromTemp(CC_temp);

            outputs = makeImageMaskedChangeData(bwfinal, bbox, necessaryData);
          
           
            outputs.spheroids_order = orderSpheroids(CC_temp, necessaryData.checkWater, necessaryData.bw_nowater);
            outputs.bwfinal = bwfinal;
            outputs.boundingBox = bbox;
            outputs.centers = centers;
            outputs.coloredLabels = label2rgb(CClabel, 'hsv', 'k', 'shuffle')*0.9;
            
end



function outputs = makeImageMaskedChangeData(bwfinal, bbox, necessaryData)
            outputs = struct;
            %Find perimeter
            bw_perim = double(bwperim(bwfinal));
            %Make visible (dilate)
            bw_perim_dil = imdilate(bw_perim,ones(3,3));
            bulkI = necessaryData.bulkI;
           
            
            
            displayI = bulkI;
            %Normalize raw image to min and max
            displayI = double(displayI);
            displayI = (displayI - min(displayI(:)))./(max(displayI(:)) - min(displayI(:)));

            imgMasked = imoverlay(displayI,bw_perim_dil, [0 0 1]);              

            [imgMasked, nonCirc, circularBoolean] = findCircularity(necessaryData, imgMasked);
            

            [imgMasked, contaminant, nonContaminantBoolean] = findContaminant(necessaryData, imgMasked, bbox);

           

            outputs.imgMasked = imgMasked;
            outputs.nonCircularCC = nonCirc;
            outputs.circularBoolean = circularBoolean;
            outputs.contaminantCC = contaminant;
            outputs.nonContaminantBoolean = nonContaminantBoolean; 


            
            
            

           

end



function [imgMaskedNew, contaminant, nonContaminantBoolean] = findContaminant(necessaryData, imgMaskedOrig, boundingBox)
            
            contaminantIntensityScale = necessaryData.contaminantIntensity;
            cc = necessaryData.cc;
            bulkI = necessaryData.bulkI;
            
            CClabel = labelmatrix(cc);
            
%             
            
            
                        
            nonContaminantBoolean = ones(1, length(cc.PixelIdxList));
                %initial struct to allow for other functions to
            %operate on it
%                 contaminant = struct;
%                 contaminant.PixelIdxList = cell(1, length(cc.PixelIdxList));

                
            
            %if contaminantIntensity ~= 0

                %Copy original cc
            contaminant = cc;
            if contaminantIntensityScale ~= 0
                
                image = bulkI;
                contaminantIntensity = contaminantIntensityScale*max(image(:));
                nonEmptySpher = 0;
               
                

              
                
                for sphere = 1:length(cc.PixelIdxList)
                    if ~isempty(cc.PixelIdxList{sphere})
                        nonEmptySpher = nonEmptySpher + 1;
                        bbox = boundingBox{nonEmptySpher};
                        imageNew = imcrop(image, bbox);
                        CCNew = imcrop(CClabel, bbox);
                        
                        wholeSpheroidPixels = double(imageNew(CCNew == sphere));
                        medianTotalPixels = median(wholeSpheroidPixels);
                        
                        
                       
                        %compare to threshold
                        if medianTotalPixels > contaminantIntensity
                            contaminant.PixelIdxList{sphere} = [];
    
                        else
                            nonContaminantBoolean(sphere) = 0;
                        end
                    end
                end
            else
                for sphere = 1:length(cc.PixelIdxList)
                    contaminant.PixelIdxList{sphere} = [];
                    nonContaminantBoolean(sphere) = 1;
                end
                
                
                
            end

                
                contaminantCCLabel = labelmatrix(contaminant);

                contaminantPerim = bwperim(contaminantCCLabel ~= 0);
                contaminant_perim_dil = imdilate(contaminantPerim, ones(3,3));
                imgMaskedNew = imoverlay(imgMaskedOrig,contaminant_perim_dil,[1 0 1]);

                
end



function [imgMaskedNew, nonCirc, circularBoolean, circ] = findCircularity(necessaryData, imgMaskedOrig)
            
            minimumCircularity = necessaryData.minimumCircularity;
            cc = necessaryData.cc;
            
                   
            circularBoolean = ones(1, length(cc.PixelIdxList));
            %initial struct to allow for other functions to
 
                

            %if minimumCircularity ~= 0
            nonCirc = cc;
            if minimumCircularity ~=0
                CC_temp = cc;
                CC_temp.PixelIdxList = cc.PixelIdxList(~cellfun('isempty',cc.PixelIdxList));
                CC_temp.NumObjects = length(CC_temp.PixelIdxList);
                %Highlight non-circular spheroids
                
                props = regionprops("table", CC_temp, ["Circularity" "Perimeter"]);
                c = props.Circularity;
                p = props.Perimeter;
                r_eq = p/(2*pi) + 0.5;
                correction = (1 - (0.5./r_eq)).^2;
                props.Circularity = min(c .* correction, 1);
                circ = cat(1,props.Circularity);           
                
                nonEmptySpher = 0;
                for sphere = 1:length(cc.PixelIdxList)
                    if ~isempty(cc.PixelIdxList{sphere})
                        nonEmptySpher = nonEmptySpher + 1;
                        if circ(nonEmptySpher) >= minimumCircularity
                            nonCirc.PixelIdxList{sphere} = [];
                            circularBoolean(sphere) = 1;
    
                        else
                            circularBoolean(sphere) = 0;
                        end
                    end
                end
            else
                for sphere = 1:length(cc.PixelIdxList)
                    nonCirc.PixelIdxList{sphere} = [];
                    circularBoolean(sphere) = 1;
                end
                ccTemp = findCCTemp(cc);
                circ = ones(length(ccTemp.PixelIdxList));
            end

                
                nonCircCCLabel = labelmatrix(nonCirc);

                nonCircPerim = bwperim(nonCircCCLabel ~= 0);
                nonCirc_perim_dil = imdilate(nonCircPerim, ones(3,3));
                imgMaskedNew = imoverlay(imgMaskedOrig,nonCirc_perim_dil,[1 1 0]);
                               
%             elseif any(~(circularBoolean))
%                 
%                 nonCircCCLabel = labelmatrix(nonCirc);
% 
%                 nonCircPerim = bwperim(nonCircCCLabel ~= 0);
%                 nonCirc_perim_dil = imdilate(nonCircPerim, ones(3,3));
%                 imgMaskedNew = imoverlay(imgMaskedOrig,nonCirc_perim_dil,[1 1 0]);
%             else
%                 imgMaskedNew = imgMaskedOrig;
%             end
            
end






function CC_temp = findCCTemp(cc)
                       
            CC_temp = cc;                     
          
            CC_temp.PixelIdxList = cc.PixelIdxList(~cellfun('isempty',cc.PixelIdxList));
            CC_temp.NumObjects = length(CC_temp.PixelIdxList);
            
        end

        %Finds bounding box from CC_temp (bounding box won't work if CC is
        %empty)
        function boundingBox = findBoundingBoxFromTemp(CC_temp)
             
            %Create bounding box matrix
            statsBW_temp = regionprops(CC_temp, 'BoundingBox');
            boundingBox = {statsBW_temp.BoundingBox};
            
        end



        function centers = findCentersFromCC(cc)
            statsBW = regionprops(cc, 'Centroid');
            centersOrig = [statsBW.Centroid];
            centers = reshape(centersOrig, 2, length(statsBW))';
            
        end

        function spheroids_order = orderSpheroids(CC_temp, checkWater, bw_nowater)            
            
            if checkWater == 1 
                %Load pre-watershed
                
                CC_nowater = bwconncomp(bw_nowater);
                
                CClabel_nowater = labelmatrix(CC_nowater);
                %Preallocate
                cluster = cell(1,CC_nowater.NumObjects+1);
                %Determine spheroid numbering
                 
                for k = 1:CC_temp.NumObjects
                    %original spheroid cluster
                   
                    location = CClabel_nowater(CC_temp.PixelIdxList{k});
                    spheroid_orig = round(mode(location(location~=0)));
                    if spheroid_orig ~= 0
                    %store cluster index
                        a = [cluster{spheroid_orig} k];
                        cluster{spheroid_orig} = a;
                    else
                        a = [cluster{end} k];
                        cluster{end} = a;
                    end
                end


                
                %Matricize and store
                spheroids_order = cell2mat(cluster);  
            else
                spheroids_order = 1:CC_temp.NumObjects;
            end
                
                
        end



        function splitBw = watershedSplit(bwToSplit, blurring, numIterations, sizeThresh)
            bw_holesTemp = bwToSplit;
            for iteration = 1:numIterations
                
                %Find distance transform and smooth
                dist = bwdist(imcomplement(double(bw_holesTemp)));
                
                if blurring > 0
                    smooth = imgaussfilt(dist,blurring);
                else
                    smooth = dist;
                end
                %smooth = dist;
                
                dist2 = -smooth;
                dist2(~bw_holesTemp) = Inf;

                %Watershed based on distance transform, splitting organoids
                L = watershed(dist2);
                bw_holesTemp(L==0) = 0;
                %Remove noise based on size threshold
                bwlil = bwareaopen(bw_holesTemp,sizeThresh);
                %Clear border
                bw_holesTemp = imclearborder(bwlil);
                bw_holesTemp = logical(bw_holesTemp);
            end
            splitBw = bw_holesTemp;
        end


        function mag = getImageMag(origI, imageReconstructionRadius, useReconstructedForGradient,  useInverted)
            if useReconstructedForGradient == 1
                se2 = strel('disk', round(imageReconstructionRadius/2));
                Iobrcbr2 = prepareImage(origI, se2, useInverted);
                [mag, ~] = imgradient(Iobrcbr2);
                    %[mag, ~] = imgradient(imgaussfilt(origI,.5));

            else
                [mag, ~] = imgradient(origI);
            end
        end


        function background  = getBackground(Iobrcbr)
            %background = imdilate(Iobrcbr, strel('Disk', round(size(Iobrcbr, 1)/20)));
            background = imfilter(Iobrcbr,fspecial('average',[500 500]),'replicate');
            
        end

        function bw_holesNew = gradientOnlyEdgeCorrect(bw_holes, edgeIntensity, Iobrcbr, bulkI, imageReconstructionRadius, useReconstructedForGradient, useInverted,sizeThresh)
            
            

            mag = getImageMag(bulkI, imageReconstructionRadius, useReconstructedForGradient,  useInverted);
            

            background = getBackground(Iobrcbr);
            mag = mag > edgeIntensity.*background;
            
            %mag = bwmorph(mag, 'spur', inf);
            bw_holes_temp = mag + bw_holes;
            
            %Remove extra area around perimeter added by
            %magnitude
            
            bwsmall2 = bwareaopen(bw_holes_temp, sizeThresh);
           
            %Smooth
            bwsmooth2 = imclose(bwsmall2,ones(3));

            %fill holes
            bw_holes_temp = bwfill(bwsmooth2,'holes');
            bw_holes_temp = imerode(bw_holes_temp, strel('disk', 3));           
            bw_holes_temp = bwareaopen(bw_holes_temp, sizeThresh);
            bw_holesNew = bw_holes+bw_holes_temp;
            
            
            
        end

        function bw_holesNew = gradientPreserveBoundaryEdgeCorrect(bw_holes, edgeIntensity, Iobrcbr, bulkI, imageReconstructionRadius, useReconstructedForGradient, useInverted,sizeThresh)
            %mark background based on segmentation
                    
                D = bwdist(bw_holes);
                DL = watershed(D);
                bgm = DL == 0;
                
                
                                  

                mag = getImageMag(bulkI,imageReconstructionRadius, useReconstructedForGradient,  useInverted);
                
                
                

                background = getBackground(Iobrcbr);

                
                %Threshold image gradient based on intensity of
                %background at that point
                mag = mag > edgeIntensity*background;
                mag = bwmorph(mag, 'spur', inf);
               
                %Add thresholded image gradient to original bwfinal
                combine = mag+bw_holes;

                %Zero-out combined bw at background marker to
                %eliminate combining of spheroids
                combine(imdilate(bgm, ones(3,3))) = 0;
                
                
                bwsmall2 = bwareaopen(combine, sizeThresh);
                
                
                bwsmooth2 = imclose(bwsmall2, ones(3));             
                bw_holes2 = imfill(bwsmooth2, 'holes');
                

               
                
                bw_holes_temp = imerode(bw_holes2, strel('disk', 3));
                bw_holes_temp = imclearborder(bw_holes_temp);
                bw_holes_temp = bwareaopen(bw_holes_temp, sizeThresh);
               

                %Remove any objects from bwfinalTemp which do not
                %appear in original bwfinal
                pixels = find(bw_holes);
                newCC = bwconncomp(bw_holes_temp);
                
                

                for spheroid = 1:length(newCC.PixelIdxList)
                    if ~any(ismember(pixels, newCC.PixelIdxList{spheroid}))
                        newCC.PixelIdxList{spheroid} = [];
                    end
                end
                
                
                finalCCLabel = labelmatrix(newCC);
                bw_holes_temp = logical(finalCCLabel);
               
                %Remove extra area around perimeter added by
                %magnitude
                %bw_holes_temp = logical(bw_holes_temp);
                
                bw_holesNew = bw_holes_temp + bw_holes;
                

                
        end

        function bw_holesNew = watershedEdgeCorrect(bw_holes, bulkI, imageReconstructionRadius, useReconstructedForGradient, useInverted,sizeThresh)
            
            D = bwdist(bw_holes);
            DL = watershed(D);
            bgm = DL == 0;         
            
            
            bw_holes_temp = bw_holes;
            bw_holes_temp(DL == 0) = 0;
            
            
            
            
            se2 = strel(ones(3,3));                        
            fgm = imerode(bw_holes_temp,se2);
            
       
            fgm = bwareaopen(fgm,sizeThresh/2);

            mag = getImageMag(bulkI, imageReconstructionRadius, useReconstructedForGradient,  useInverted);

            %labels = imdilate(bgm,ones(3,3)) + 3*fgm;

            %tempMag = double(mag);
            %tempMag = (tempMag - min(tempMag(:)))./(max(tempMag(:)) - min(tempMag(:)));

             
            
            
            
            
            mag = imgaussfilt(mag,1);
            mag2 = imimposemin(mag,  fgm | bgm );
            
            L = watershed(mag2);
            
            mask = imfill(L==0, 'holes');
            combine = mask+bw_holes_temp;

            bwsmall2 = bwareaopen(combine, sizeThresh);                       
            %bw_holesNew = imclearborder(bwsmall2);
            bw_holesNew = bwsmall2;

        end