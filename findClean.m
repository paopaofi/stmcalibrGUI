function findClean()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Im=sxmopen("/Users/bubblefish/Documents/CDMPhy/project ab initio/REPORT1/Au/Au/Bi2201#19_460.sxm");
    %Im=sxmopen("/Users/bubblefish/Documents/CDMPhy/project ab initio/REPORT1/Au/Au/Au_20200719_01_4.2K_011.sxm");
    testingIm = Im{1};
    Rscan = imref2d(size(testingIm.data));
    Rscan.XWorldLimits = [0, testingIm.width];
    Rscan.YWorldLimits = [0, testingIm.height];
    ImOrigin = testingIm.data;
    ImFiltered = Im_Flatten_XY(ImOrigin);
    mask = FindSteps(ImFiltered);
    ImFs = Im_Flatten_XY(ImOrigin,mask);
    ImFs = (ImFs - min(ImFs(:)))/(max(ImFs(:)) - min(ImFs(:)));
    ImGauF1 = imgaussfilt(ImFs,5.2);
    ImGauF2 = imgaussfilt(ImFs,3.1);
    ImGauF3 = imgaussfilt(ImFs,1);
    med = median(ImFs(:));
    
    figure();
    
    subplot(3,3,1)
    imshow(ImFs,Rscan,[]);
    
    % 1-d Image histogram analysis part
    % load data set
    binWidth = 0.003333;
    binCtrs = 0:binWidth:1;
    PixeLum = ImFs(:);
    nu = length(PixeLum);
    counts = hist(PixeLum,binCtrs);
    prob = counts/(nu*binWidth);
    
    % plot normalized histogram
    % plot(binCtrs, prob, 'b-', 'linewidth', 1.5);
    subplot(3,3,2);
    plot(binCtrs, prob, 'ko');
    hold on;
    
    % find a serise of multigaussian models
    ComNum = 1:12;
    ComNum = 3*ComNum;
    AIC = zeros(1,12);
    BIC = zeros(1,12);
    ConVrg = zeros(1,12);
    for i = 1:12
        obj{i} = gmdistribution.fit(PixeLum,i);
        AIC(i) = obj{i}.AIC;
        BIC(i) = obj{i}.BIC;
        ConVrg(i) = obj{i}.Converged;
    end
    CAIC = BIC + ComNum; % calculate the CAIC information criterion
    % find the optimum number of gaussian components
    [~,numComponent_C] = min(CAIC);
    
    % generate corresponding pdfs
    paramEsts_C = obj{numComponent_C};
    MU_C = paramEsts_C.mu;
    SIGMA_C = paramEsts_C.Sigma;
    PPp_C = paramEsts_C.PComponents;
    OptimumFit_C = gmdistribution(MU_C,SIGMA_C,PPp_C);
    
    % plot best models
    xgridss = transpose(linspace(0,1,301));
    plot(xgridss,pdf(OptimumFit_C,xgridss),'-','linewidth',3);
    for i = 1:numComponent_C
        plot(xgridss,PPp_C(i)*pdf(gmdistribution(MU_C(i),SIGMA_C(i),1),xgridss),'-','linewidth',1)
    end
    
    % determined the 3 sigma threshold luminosity for the highest terrace
    ParaMap = find(PPp_C > 0.05);
    [maxMU, indexMU] = max(MU_C(ParaMap));
    maxSIG = sqrt(SIGMA_C(ParaMap(indexMU)));
    Threshold = maxMU + 3*maxSIG;
    plot([Threshold Threshold],[0 max(prob)],'r');
    hold off;
    
    % Screen the original image with Luminosity higher than the Threshold value
    ImFree = ImFs;
    thrsMask = (ImFs > Threshold);
    ImFree(thrsMask(:)) = med;
    % re-flatten and re-normalization
    % ImEdges = FindSteps(ImFree);
    % ImFree = Im_Flatten_XY(ImFree,ImEdges);
    % ImFree = (ImFree - min(ImFree(:)))/(max(ImFree(:)) - min(ImFree(:)));
    ImEdges = FindSteps(ImFree);
    
    % plot screened image
    subplot(3,3,3);
    imshow(ImFree, Rscan, []);
    
    % plot the edge features
    subplot(3,3,4);
    imshow(ImEdges,Rscan,[]);
    
    % find on-edge dirty spots using branchpoint detection and thicken
    % operation
    DirtyonEdge = bwmorph(ImEdges,'branchpoints');
    DirtyonEdge = bwmorph(DirtyonEdge,'thicken');
    DirtyonEdge = DirtyonEdge|bwmorph(bwmorph(FindSteps(ImFs),'branchpoints'),'thicken');
    subplot(3,3,5);
    imshow(DirtyonEdge,Rscan,[]);
    
    
    % small dirty spots detection using a combination of edge and imfill
    bu = edge(ImFs, 'Canny',0.02);
    bu = bwmorph(bu,'hbreak');
    buGau1 = edge(ImGauF1,'Canny',0.08);
    buGau2 = edge(ImGauF2,'Canny',0.08);
    buGau3 = edge(ImGauF3,'Canny',0.02);
    DirtySpots = imfill(bu,'holes') - bu;
    SpotsGau = imfill(buGau1,'holes')-buGau1;
    SpotsGau = SpotsGau|(imfill(buGau2,'holes')-buGau2);
    SpotsGau = SpotsGau|(imfill(buGau3,'holes')-buGau3);
    CombinedSpots = DirtySpots|SpotsGau;
    subplot(3,3,6);
    imshow(bu, Rscan, [], 'Border', 'tight');
    subplot(3,3,7);
    imshow(CombinedSpots,Rscan,[]);
    
    % display all the features detected
    subplot(3,3,8)
    dirtymask = CombinedSpots|DirtyonEdge;
    dirtymask = dirtymask|thrsMask;
    allmask = dirtymask|ImEdges;
    imshow(allmask,Rscan,[]);
    
    % display the calculated clean reigons 
    subplot(3,3,9)
    xnumber = round(testingIm.width/5e-9);
    ynumber = round(testingIm.height/5e-9);
    [yres,xres] = size(ImFs);
    xpixelNum = xres/xnumber;
    ypixelNum = yres/ynumber;
    imshow(ImFs,Rscan,[]);
    
    
    for i = 1:xnumber
        for j = 1:ynumber
            location = [(i-0.5)*5e-9, (j-0.5)*5e-9, 5e-09, 5e-09];
            EdgePixels = sum(sum(ImEdges((j-1)*ypixelNum+2:j*ypixelNum-1, (i-1)*xpixelNum+2:i*xpixelNum-1)));
            DirtyPixels = sum(sum(dirtymask((j-1)*ypixelNum+2:j*ypixelNum-1, (i-1)*xpixelNum+2:i*xpixelNum-1)));            
            if DirtyPixels < 1
                rg = rectangle(gca,'Position',[location(1)-0.5*location(3),location(2)-0.5*location(4),location(3),location(4)],'EdgeColor','g');
                rg.LineWidth = 0.3;
                if EdgePixels > ((xpixelNum + ypixelNum)/4) 
                    rg.EdgeColor = 'r';
                    rg.LineWidth = 2.0;
                end
            end
        end
    end
    
    
end

