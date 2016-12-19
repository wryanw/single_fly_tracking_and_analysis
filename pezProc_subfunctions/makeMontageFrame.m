function [imData,outcome] = makeMontageFrame(vidName,savePath,Style)
imData = [];
outcome = 'pass';
dwnSampl = 1;%can be reduced to slightly improve speed at the cost of quality

if nargin == 0
%     outcome = 'no input received';
%     return
    vidName = 'run010_pez3003_20140312_expt0008000002590057_vid0016.mp4';
    savePath = 'C:\Users\cardlab\Documents\pictest.tif';
    Style = 1;
elseif nargin == 1
    savePath = 'none';
    Style = 1;%default is top, side view, zoomed in
elseif nargin == 2
    Style = 1;
end

%%%%% computer and directory variables and information
op_sys = system_dependent('getos');
if strfind(op_sys,'Microsoft Windows 7')
    archDir = [filesep filesep 'arch' filesep 'card'];
    dm11Dir = [filesep filesep 'dm11' filesep 'cardlab'];
else
    archDir = [filesep 'Volumes' filesep 'card'];
    if ~exist(archDir,'file')
        archDir = [filesep 'Volumes' filesep 'card-1'];
    end
    dm11Dir = [filesep 'Volumes' filesep 'cardlab'];
end
if ~exist(archDir,'file')
    error('Archive access failure')
end
if ~exist(dm11Dir,'file')
    error('dm11 access failure')
end
parentDir = fullfile(archDir,'Data_pez3000');


strParts = strsplit(vidName,'_');
runDir = [strParts{1} '_' strParts{2} '_' strParts{3}];
runPath = fullfile(parentDir,strParts{3},runDir);
vidStatsPath = fullfile(runPath,[runDir '_videoStatistics.mat']);
if ~exist(vidStatsPath,'file')
    outcome = 'Video statistics does not exist';
    return
end
vidStats_loading = load(vidStatsPath);
datanameLoad = fieldnames(vidStats_loading);
vidStats = vidStats_loading.(datanameLoad{1});
vidPath = fullfile(runPath,vidName);
[~,videoID,vidExt] = fileparts(vidPath);
if ~exist(vidPath,'file')
    outcome = 'Video not found';
    return
end
try
    downloadOpRestrictedTest = strcmp(vidStats.download_option{videoID},'Restricted full rate');
    downloadOpCutrateTest = strcmp(vidStats.download_option{videoID},'Save Cut Rate');
    downloadOpFullrateTest = strcmp(vidStats.download_option{videoID},'Save Full Rate');
catch 
    outcome = 'Unexpected video statistics data';
    return
end
if downloadOpRestrictedTest
    vidNameSupp = [videoID '_supplement' vidExt];
    vidSuppPath = fullfile(runPath,'highSpeedSupplement',vidNameSupp);
    if ~exist(vidSuppPath,'file')
        outcome = 'Supplement does not exist';
        return
    end
    vidObj = VideoReader(vidSuppPath);
    suppFrmCt = vidObj.NumberOfFrames;
    suppTargetCt = numel(vidStats.supplement_frame_reference{videoID});
    if suppFrmCt ~= suppTargetCt
        outcome = 'Unexpected supplement frame count';
        return
    end
end
vidObj = VideoReader(vidPath);
recFrmCt = vidObj.NumberOfFrames;
if downloadOpCutrateTest || downloadOpRestrictedTest
    masterFrmRefs = vidStats.cutrate10th_frame_reference{videoID};
elseif downloadOpFullrateTest
    masterFrmRefs = (1:vidStats.frame_count(videoID));
else
    outcome = 'Unknown download option';
    return
end
recTargetCt = numel(masterFrmRefs);
if recFrmCt ~= recTargetCt
    outcome = 'Unexpected cutrate frame count';
    return
end

vidH = vidObj.Height;
vidW = vidObj.Width;
roiPos = vidStats.roi{videoID};
stagePos = vidStats.prism_base{videoID};
avgPrismL = mean([roiPos(3)-roiPos(1),roiPos(4)-roiPos(2)]);

switch Style
    case 1 %side view, zoomed in
        dim1b = mean(stagePos(:,2))*1.1;
        dim1a = dim1b-avgPrismL/2;
        if dim1a < 1, dim1a = 1; end
        dwnAa = round(linspace(dim1a,dim1b,round(dwnSampl*(dim1b-dim1a))))';
        dim2a = vidW/2-avgPrismL/2;
        if dim2a < 1, dim2a = 1; end
        dim2b = vidW/2+avgPrismL/2;
        if dim2b > vidW, dim2b = vidW; end
        dwnBa = round(linspace(dim2a,dim2b,round(dwnSampl*(dim2b-dim2a))));
        dwnAb = repmat(dwnAa,1,numel(dwnBa));
        dwnBb = repmat(dwnBa,numel(dwnAa),1);
        dwnC = dwnAb+(dwnBb-1).*vidH;
        dwnS = size(dwnC);
    case 2
        dim1b = mean(stagePos(:,2))*1.1;
        dim1a = 1;
        dwnAa = round(linspace(dim1a,dim1b,round(dwnSampl*(dim1b-dim1a))))';
        dim2a = 1;
        dim2b = vidW;
        if dim2b > vidW, dim2b = vidW; end
        dwnBa = round(linspace(dim2a,dim2b,round(dwnSampl*(dim2b-dim2a))));
        dwnAb = repmat(dwnAa,1,numel(dwnBa));
        dwnBb = repmat(dwnBa,numel(dwnAa),1);
        dwnC = dwnAb+(dwnBb-1).*vidH;
        dwnS = size(dwnC);
    case 3
        dim1b = roiPos(4);
        dim1a = roiPos(2);
        dwnAa = round(linspace(dim1a,dim1b,round(dwnSampl*(dim1b-dim1a))))';
        dim2a = roiPos(1);
        dim2b = roiPos(3);
        if dim2b > vidW, dim2b = vidW; end
        dwnBa = round(linspace(dim2a,dim2b,round(dwnSampl*(dim2b-dim2a))));
        dwnAb = repmat(dwnAa,1,numel(dwnBa));
        dwnBb = repmat(dwnBa,numel(dwnAa),1);
        dwnC = dwnAb+(dwnBb-1).*vidH;
        dwnS = size(dwnC);
    otherwise
        outcome = 'unknown style';
        return
end
rowCt = 5;
colCt = 4;
montyRefs = round(linspace(1,recFrmCt,rowCt*colCt));
montyIm = uint8(zeros(dwnS.*[rowCt colCt]));
frmRef = cell(rowCt*colCt,1);
labelPos = cell(rowCt*colCt,1);
for iterR = 1:rowCt
    for iterC = 1:colCt
        baseRef = (iterR-1)*colCt+iterC;
        frmRead = read(vidObj,montyRefs(baseRef));
        frmRead = frmRead(:,:,1);
        colRef = ((iterC-1)*dwnS(2)+1:(iterC)*dwnS(2));
        rowRef = ((iterR-1)*dwnS(1)+1:(iterR)*dwnS(1));
        montyFrm = frmRead(dwnC);
        montyFrm(1,:) = 50;
        montyFrm(end,:) = 50;
        montyFrm(:,1) = 50;
        montyFrm(:,end) = 50;
        montyIm(rowRef,colRef) = montyFrm;
        frmRef{baseRef} = masterFrmRefs(montyRefs(baseRef));
        labelPos{baseRef} = [colRef(1)+5 rowRef(1)+30]./fliplr(dwnS.*[rowCt colCt]);
    end
end
imAdjInput = stretchlim(montyIm,[0.01 0.999]);
imData = repmat(imadjust(montyIm,imAdjInput),[1 1 3]);
for iterT = 1:rowCt*colCt
        imData = textN2im_v2(imData,{num2str(frmRef{iterT})},14,labelPos{iterT},'left',12,0);
end

if ~strcmp(savePath,'none')
    imwrite(imData,savePath)
    imData = [];
% else
%     imshow(imData)
%     imData = [];
end


% end