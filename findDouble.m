function findDouble()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Im=sxmopen("/Users/bubblefish/Documents/CDMPhy/project ab initio/REPORT1/Au/Au/Bi2201#19_460.sxm");
    %Im=sxmopen("/Users/bubblefish/Documents/CDMPhy/project ab initio/REPORT1/Au/Au/Au_20200719_01_4.2K_002.sxm");
    testingIm = Im{1};
    Rscan = imref2d(size(testingIm.data));
    Rscan.XWorldLimits = [0, testingIm.width];
    Rscan.YWorldLimits = [0, testingIm.height];
    ImOrigin = testingIm.data;
    ImFiltered = Im_Flatten_XY(ImOrigin);
    mask = FindSteps(ImFiltered);
    ImFs = Im_Flatten_XY(ImOrigin,mask);
    ImFs = (ImFs - min(ImFs(:)))/(max(ImFs(:)) - min(ImFs(:)));
    
    figure();

    scanLine = redRec(20,:);
    Dscan = abs(diff(scanLine));

    [~, num] = size(Dscan);
    scansum = sum(Dscan);

