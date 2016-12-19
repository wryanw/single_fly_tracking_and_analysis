function graphTable = makeGraphingTable
%%

if isempty(mfilename)
    clearvars
end

%%%%%% MING PAPER
excelDir = 'Z:\Ming_RubinLab';
% excelName = 'Table1_v5_WRW_MW_v2_swapped.xlsx';
% excelName = 'backing_on_pez.xlsx';
excelName = 'jumping_on_pez.xlsx';
% excelName = 'inactivation_table.xlsx';
% excelName = 'mingPaperFinalList.xlsx';

% excelDir = 'Z:\Ming_RubinLab\FlyPez_backing';
% excelName = 'DL_looming.xlsx';

%%%%%% METHODS PAPER
% excelDir = 'Z:\pez3000_methods_paper';
% excelName = 'expts_in_mings_paper.xlsx';
% excelName = 'methods_paper_activation.xlsx';
% excelName = 'fixing.xlsx';
% excelName = 'low_activation_DI.xlsx';
% excelDir = 'Z:\Denise';
% excelName = 'DI_rescreen_lines.xlsx';

% excelName = 'DL_psychophys.xlsx';
% excelName = 'OL_screen_methodsPaper.xlsx';
% excelName = 'LC4_LC6.xlsx';
% excelName = 'LC4_LC6_GF_DL.xlsx';

% excelDir = 'Z:\WryanW\flyPez_experiment_management';
% excelName = 'CX_screen.xlsx';
% excelName = 'LPLC2_exploration.xlsx';
% excelDir = 'C:\Users\williamsonw\Documents';
% excelName = 'DI.xlsx';

% excelDir = 'Z:\Martin';
% excelName = 'playingwithryan.xlsx';

excelPath = fullfile(excelDir,excelName);
makeGraphOptionsStruct(excelPath)
%%
dataTable = loadDataTable;
optionsPath = 'Z:\Data_pez3000_analyzed\WRW_graphing_variables\graphOptions.mat';
graphOptions = load(optionsPath);
graphOptions = graphOptions.graphOptions;
useManual = graphOptions.useManual;
excelPath = graphOptions.excelPath;
sheetName = graphOptions.exptSheet;
sliceOption = graphOptions.sliceOption;
excelTable = readtable(excelPath,'ReadRowNames',true,'Sheet',sheetName);
excelVarNames = excelTable.Properties.VariableNames;
if max(strcmp(excelVarNames,graphOptions.unique_label)) == 0
    error('specified unique_label not found in excel file exptIDsheet')
end
excelVarNames{strcmp(excelVarNames,graphOptions.unique_label)} = 'unique_label';
excelTable.Properties.VariableNames = excelVarNames;
exptIDlist = excelTable.Properties.RowNames;
exptIDlist = cellfun(@(x) x(4:end),exptIDlist,'uniformoutput',false);
excelTable(cellfun(@(x) numel(x),exptIDlist) ~= 16,:) = [];
exptIDlist = excelTable.Properties.RowNames;
exptCt = numel(exptIDlist);
graphIDlist = unique(excelTable.unique_label);
graphCt = numel(graphIDlist);

vidList = dataTable.Properties.RowNames;
dataExptIDs = cellfun(@(x) x(29:44),vidList,'uniformoutput',false);
for iterE = 1:exptCt
    exptID = exptIDlist{iterE}(4:end);
    labelNdx = strcmp(dataExptIDs,exptID);
    dataTable.unique_label(labelNdx) = repmat(excelTable.unique_label(['ID#' exptID]),sum(labelNdx),1);
    if sliceOption > 1
        if max(strcmp(excelTable.Properties.VariableNames,'sliceOp_azi'))
            dataTable.sliceOp(labelNdx) = repmat(excelTable.sliceOp_azi(['ID#' exptID]),sum(labelNdx),1);
        else
            dataTable.sliceOp(labelNdx) = ones(sum(labelNdx),1);
        end
    elseif sliceOption == -1
        if max(strcmp(excelTable.Properties.VariableNames,'sliceOp_lv'))
            dataTable.sliceOp(labelNdx) = repmat(excelTable.sliceOp_lv(['ID#' exptID]),sum(labelNdx),1);
        else
            dataTable.sliceOp(labelNdx) = ones(sum(labelNdx),1);
        end
    end
    if strcmp(excelTable.Experiment_Type{['ID#' exptID]},'Combo')
        dataTable.stimInfo(labelNdx) = repmat(excelTable.activation1(['ID#' exptID]),sum(labelNdx),1);
    end
end
%%
if sliceOption == 1
    if ~max(strcmp(excelTable.Properties.VariableNames,'activation2'))
        error('activation2 is not available')
    end
    actList = unique(cat(1,excelTable.activation1,excelTable.activation2,excelTable.activation3));
    actList(strcmp(actList,'None')) = [];
    if ~(numel(unique(excelTable.activation1)) > 1 || numel(actList) > 1)
        error('multiple stimuli not detected')
    end
    ids2keep = strcmp(excelTable.Experiment_Type,'Photoactivation');
    ids2keep = strcmp(excelTable.Experiment_Type,'Combo') | ids2keep;
    activationIDs = exptIDlist(ids2keep);
    activationIDs = cellfun(@(x) x(4:end),activationIDs,'uniformoutput',false);
    vidList = dataTable.Properties.RowNames;
    dataExptIDs = cellfun(@(x) x(29:44),vidList,'uniformoutput',false);
    activationNdx = cellfun(@(x) max(strcmp(activationIDs,x)),dataExptIDs);
    try
        actKey = readtable(excelPath,'ReadRowNames',false,'Sheet','activation_key');
    catch
        try
            actList = cat(1,{'None'},actList(:));
            newList = cat(1,{'unknown'},cell(numel(actList)-1,1));
            writetable(cell2table([actList newList],'VariableNames',{'oldName','newName'}),excelPath,...
                'Sheet','activation_key','WriteVariableNames',true)
        catch
            disp('Please close excel file and rerun')
        end
        disp({'Please complete the activation key using ''unknown'' as newName',...
            'for protocols to omit and then rerun code'}')
        return
    end
    actKey = table2cell(actKey);
    actKey = actKey(:,1:2);
    for iterR = 1:size(actKey,1)
        actIndex = cellfun(@(x) ~isempty(strfind(x,actKey{iterR,1})),dataTable.stimInfo);
        dataTable.stimInfo(actIndex) = repmat(actKey(iterR,2),sum(actIndex),1);
    end
    dataTable.unique_label(activationNdx) = cellfun(@(x,y) cat(2,x,'_',y),dataTable.unique_label(activationNdx),dataTable.stimInfo(activationNdx),'uniformoutput',false);
    dataTable.unique_label(~activationNdx) = cellfun(@(x,y) cat(2,x,'_noActivation'),dataTable.unique_label(~activationNdx),'uniformoutput',false);
    dataTable(strcmp(dataTable.stimInfo,'unknown'),:) = [];
    graphIDlist = unique(dataTable.unique_label);
    
elseif sliceOption > 1
    aziGrpCt = sliceOption;% 5 or 7 is recommended ... 13 works too
%     if aziGrpCt/2 == round(aziGrpCt/2)
%         aziBreaks = linspace(0,180,aziGrpCt*2+1);
%         aziNames = aziBreaks(2:2:aziGrpCt*2+1);
%         aziBreaks = aziBreaks(1:2:aziGrpCt*2+1);
%     else
        aziBreaks = linspace(0,180,aziGrpCt*2-1);
        aziNames = aziBreaks(1:2:aziGrpCt*2-1);
        aziBreaks = [aziBreaks(1) aziBreaks(1,2:2:aziGrpCt*2-1) aziBreaks(end)];
        if ~isempty(strfind(graphOptions.sheetName,'fieldMap'))
            overlap = median(diff(aziNames))/2;
        else
            overlap = 0;
        end
%     end
    overCell = cell(graphCt,aziGrpCt);
    for iterG = 1:graphCt
        sliceBool = dataTable.sliceOp;
        stimInit = dataTable.zeroFly_StimAtStimStart;
        graphID = graphIDlist{iterG};
        groupLogical = strcmp(dataTable.unique_label,graphID);
        
        for iterA = 1:aziGrpCt
            aziNdcs = stimInit > aziBreaks(iterA) & stimInit <= aziBreaks(iterA+1) & groupLogical & sliceBool;
            graphIDcell = repmat({[graphID '_azi' num2str(aziNames(iterA))]},sum(aziNdcs),1);
            dataTable.unique_label(aziNdcs) = graphIDcell;
            aziNdcsOver = stimInit > aziBreaks(iterA)-overlap & stimInit <= aziBreaks(iterA+1)+overlap & groupLogical & sliceBool;
            aziNdcsOver = aziNdcsOver & ~aziNdcs;
            dataTableOver = dataTable(aziNdcsOver,:);
            dataTableOver.Properties.RowNames = cellfun(@(x) cat(2,x,'_overlap'),dataTableOver.Properties.RowNames,'uniformoutput',false);
            graphIDcell = repmat({[graphID '_azi' num2str(aziNames(iterA))]},sum(aziNdcsOver),1);
            dataTableOver.unique_label = graphIDcell;
            overCell{iterG,iterA} = dataTableOver;
        end
        groupLogical = strcmp(dataTable.unique_label,graphID);
        dataTable(groupLogical & sliceBool,:) = [];
    end
    dataTable = [dataTable;cat(1,overCell{:})];
    graphIDlist = unique(dataTable.unique_label);
elseif sliceOption == -1
    fot = dataTable.manualFot;
    fowu_fowd_folp = cell2mat(dataTable.manual_wingup_wingdwn_legpush);
    fowu = fowu_fowd_folp(:,1);
    shortThresh = 41;
    shortLogical = (fot-fowu) <= shortThresh;%41
    longLogical = (fot-fowu) > shortThresh;% 30 200
    longLogical(isnan(longLogical)) = false;
    shortLogical(isnan(shortLogical)) = false;
    if max(strcmp(dataTable.Properties.VariableNames,'sliceOp'))
        sliceBool = dataTable.sliceOp;
    else
        sliceBool = true(size(dataTable,1),1);
    end
    shortLogical = shortLogical & sliceBool;
    longLogical = longLogical & sliceBool;
    dataTable.unique_label(shortLogical) = cellfun(@(x) cat(2,x,'_short'),dataTable.unique_label(shortLogical),'uniformoutput',false);
    dataTable.unique_label(longLogical) = cellfun(@(x) cat(2,x,'_long'),dataTable.unique_label(longLogical),'uniformoutput',false);
    dataTable((~(shortLogical | longLogical) & sliceBool),:) = [];
    graphIDlist = unique(dataTable.unique_label);
    
end
graphIDlist(cellfun(@(x) ~isempty(strfind(x,'omit')),graphIDlist)) = [];
graphCt = numel(graphIDlist);

% table variables that might not save to excel
videoList = cell(graphCt,1);
exptIDs = cell(graphCt,1);
jumpTest = cell(graphCt,1);
shortTest = cell(graphCt,1);
moveTest = cell(graphCt,1);
nonMovers = cell(graphCt,1);
earlyMovers = cell(graphCt,1);
fowu_fowd_folp = cell(graphCt,1);
jumpEnd = cell(graphCt,1);
trkEnd = cell(graphCt,1);
maxVelInfo = cell(graphCt,1);
stimInfo = cell(graphCt,1);
stimStart = cell(graphCt,1);
zeroFly_XYZmm_Tdeg_fac1000 = cell(graphCt,1);
dist_speed_accel_fac100 = cell(graphCt,1);
zeroFly_Trajectory = cell(graphCt,1);
zeroFly_Jump = cell(graphCt,1);
zeroFly_Departure = cell(graphCt,1);
zeroFly_StimAtStimStart = cell(graphCt,1);
zeroFly_StimAtJump = cell(graphCt,1);
zeroFly_StimAtFrmOne = cell(graphCt,1);
rawFly_XYZpix_Trad_frmOne = cell(graphCt,1);
distFromCenter_pix = cell(graphCt,1);
relMotion_FB_LR_Tdeg_fac100 = cell(graphCt,1);
relPosition_FB_LR_Tdeg_fac100 = cell(graphCt,1);
lvList = cell(graphCt,1);
flyLength_mm = cell(graphCt,1);
groupID = cell(graphCt,1);
plotID = cell(graphCt,1);

% table variable that can save to excel
dataCount = zeros(graphCt,1);
vidCount = zeros(graphCt,1);
plotCount = zeros(graphCt,1);
jumpCount = zeros(graphCt,1);
manCount = zeros(graphCt,1);
recordRate = zeros(graphCt,1);
stimDur = zeros(graphCt,1);
elevation = zeros(graphCt,1);
azimuth = zeros(graphCt,1);
startSize = zeros(graphCt,1);
stopSize = zeros(graphCt,1);
lv = zeros(graphCt,1);
contrast = cell(graphCt,1);
width = zeros(graphCt,1);
frequency = zeros(graphCt,1);
genotype = cell(graphCt,1);
stimType = cell(graphCt,1);
visStimType = cell(graphCt,1);
order = zeros(graphCt,1);
repeat = zeros(graphCt,1);

graphTable = table(elevation,azimuth,startSize,stopSize,lv,contrast,width,frequency,genotype,...
    vidCount,dataCount,plotCount,jumpCount,manCount,recordRate,stimDur,videoList,exptIDs,...
    jumpTest,shortTest,moveTest,nonMovers,earlyMovers,fowu_fowd_folp,jumpEnd,trkEnd,stimInfo,...
    maxVelInfo,stimStart,zeroFly_XYZmm_Tdeg_fac1000,zeroFly_Trajectory,...
    zeroFly_Jump,zeroFly_Departure,zeroFly_StimAtStimStart,zeroFly_StimAtJump,...
    zeroFly_StimAtFrmOne,rawFly_XYZpix_Trad_frmOne,dist_speed_accel_fac100,distFromCenter_pix,...
    relMotion_FB_LR_Tdeg_fac100,relPosition_FB_LR_Tdeg_fac100,lvList,groupID,...
    plotID,order,repeat,stimType,visStimType,flyLength_mm,'RowNames',graphIDlist);

for iterG = 1:graphCt
    graphID = graphIDlist{iterG};
    groupLogical = strcmp(dataTable.unique_label,graphID);
    exptTable = dataTable(groupLogical,:);
    if abs(sliceOption) > 0
        if sliceOption == 1
            strSplitRefs = strfind(graphID,'_');
            compareID = graphID(1:strSplitRefs(end)-1);
        elseif sliceOption == -1
            compareID = regexprep(graphID,'_long','');
            compareID = regexprep(compareID,'_short','');
        else
            strSplitRefs = strfind(graphID,'_azi');
            if isempty(strSplitRefs)
                compareID = graphID;
            else
                compareID = graphID(1:strSplitRefs(end)-1);
            end
        end
        smlExcelTable = excelTable(strcmp(excelTable.unique_label,compareID),:);
    else
        smlExcelTable = excelTable(strcmp(excelTable.unique_label,graphID),:);
    end
    stimTypeString = smlExcelTable.Experiment_Type{1};
    graphTable.visStimType(graphID) = smlExcelTable.visStimType(1);
    if strcmp('Visual_stimulation',stimTypeString) || strcmp('Combo',stimTypeString)
        querryName = {'azimuth','elevation','startSize','stopSize','lv','contrast','width','frequency'};
        members = ismember(querryName,smlExcelTable.Properties.VariableNames);
        for iterM = 1:numel(members)
            querryVal = smlExcelTable.(querryName{iterM})(1);
            if members(iterM) == 0, continue, end
            if strcmp(querryName{iterM},'startSize') && iscell(smlExcelTable.(querryName{iterM})(1))
                querryVal = str2double(regexprep(querryVal{1},'pt','\.'));
            end
            graphTable.(querryName{iterM})(graphID) = querryVal;
        end
        if sliceOption > 1 && exptTable.sliceOp(1) == 1
            exptAzi = graphID;
            exptAzi(1:strfind(exptAzi,'azi')+2) = [];
            exptAzi = str2double(exptAzi);
            graphTable.azimuth(graphID) = exptAzi;
        end
    end
    graphTable.recordRate(graphID) = smlExcelTable.Record_Rate(1);
    if max(strcmp(smlExcelTable.Properties.VariableNames,'genotype'))
        graphTable.genotype(graphID) = smlExcelTable.genotype(1);
    end
    graphTable.vidCount(graphID) = sum(smlExcelTable.Total_Videos);
    graphTable.stimType(graphID) = {stimTypeString};
    smlListIDs = smlExcelTable.Properties.RowNames;
    graphTable.exptIDs(graphID) = {smlListIDs};
    stimDur = max(exptTable.stimDuration(exptTable.stimDuration < Inf));
    if isinf(stimDur) || isempty(stimDur) || isnan(stimDur)
        error('stim dur error')
    end
    graphTable.stimDur(graphID) = stimDur;
    
    noStimNdcs = cellfun(@(x) ~isempty(strfind(x,'whiteonwhite')),exptTable.stimInfo);
    stimStart = exptTable.stimStart;
    if ~isempty(strfind(graphID,'noStim'))
        stimStart(noStimNdcs) = exptTable.trkEnd(noStimNdcs);
        stimInit = exptTable.zeroFly_StimAtFrmOne;
    else
        stimInit = exptTable.zeroFly_StimAtStimStart;
    end
    if strcmp('Visual_stimulation',stimTypeString)
        %%%%% Was the stimulus within acceptable boundaries??
        exptAzi = zeros(size(exptTable,1),1);
        for iterLV = 1:size(exptTable,1)
            stimInfo = exptTable.stimInfo{iterLV};
            stimInfo(1:strfind(stimInfo,'azi')+2) = [];
            stimInfo(strfind(stimInfo,'_'):end) = [];
            exptAzi(iterLV) = str2double(stimInfo);
        end
    else
        exptAzi = 0;
    end
    if isempty(exptTable)
        continue
    end
    stimTol = graphOptions.azimuth_tolerance;
    stimTest = abs(stimInit-exptAzi) < stimTol;
    if isinf(stimTol)
        stimTest = true(size(stimTest));
    end
    notCounting = ~stimTest;
    
    sizeTest = cellfun(@(x) size(x,2) ~= 3,exptTable.manual_wingup_wingdwn_legpush);
    exptTable.manual_wingup_wingdwn_legpush(sizeTest) = repmat({NaN(1,3)},sum(sizeTest),1);
    if useManual == 1
        jumpTest = exptTable.manualJumpTest;
        fot = exptTable.manualFot;
        fowu_fowd_folp = cell2mat(exptTable.manual_wingup_wingdwn_legpush);
        fowu = fowu_fowd_folp(:,1);
        shortLogical = (fot-fowu) <= 41;
    elseif useManual == 2
        fot = exptTable.autoFot;
        fowu_fowd_folp = NaN(size(fot,1),3);
        shortLogical = NaN(size(fot,1));
        jumpTest = exptTable.autoJumpTest;
    elseif useManual == 3
        jumpTest = exptTable.manualJumpTest;
        fot = exptTable.manualFot;
        fot(isnan(fot)) = exptTable.autoFot(isnan(fot));
        jumpTest(isnan(jumpTest)) = exptTable.autoJumpTest(isnan(jumpTest));
        fowu_fowd_folp = cell2mat(exptTable.manual_wingup_wingdwn_legpush);
        fowu = fowu_fowd_folp(:,1);
        shortLogical = (fot-fowu) <= 41;
    end
    
    notCounting(isnan(jumpTest)) = true;
    jumpTest(isnan(jumpTest)) = 0;
    jumpTest = logical(jumpTest);
    nonJumpers = ~jumpTest;
    nonJumpers(notCounting) = false;
    jumpTest(notCounting) = false;
    
    frms_per_ms = double(smlExcelTable.Record_Rate(1)/1000);
    ms2frm = @(x) round(x*frms_per_ms);
    moveTest = false(size(exptTable,1),1);
    nonMovers = false(size(exptTable,1),1);
    earlyMovers = false(size(exptTable,1),1);
    maxVelInfo = NaN(size(exptTable,1),5);
    moveThresh = 0.025;% 25 mm/s set by referring to robie et al 2010, Fig.3 JExpBiol
    % and based on trough seen in nonjumping max vel plot 'maxMotion'
    stillThresh = 0.005;
    distCapThresh = 1; % mm to travel before cutoff
    if strcmp('Photoactivation',stimTypeString) || ~strcmp(smlExcelTable.visStimType{1},'loom')
        mvntThreshFrms = 1;
    else
        lv = graphTable.lv(graphID);
        degA = graphTable.startSize(graphID);
        degB = 20;
        mvntThreshMS = lv/tan(deg2rad(degA/2))-lv/tan(deg2rad(degB/2));%duration from stim start to thresh size
        mvntThreshFrms = ms2frm(mvntThreshMS+20);% adding neural delay 20 ms
        if mvntThreshFrms > graphTable.stimDur(graphID)*0.75
            mvntThreshFrms = round(graphTable.stimDur(graphID)*0.75);
        end
        if mvntThreshFrms < 1
            mvntThreshFrms = 1;
        end
    end
    winLeg = ms2frm(5);
    for iterV = 1:size(exptTable,1)
        mvntArray = double(exptTable.dist_speed_accel_fac100{iterV})/100;
        mvntVec = mvntArray(stimStart(iterV):end,2);
        if numel(mvntVec) < winLeg*2+1, continue, end
        mvntVec = smooth(mvntVec,winLeg*2+1);
        mvntVec(1:winLeg) = NaN;
        mvntVec(end-winLeg:end) = NaN;
        if mvntThreshFrms > numel(mvntVec), continue, end
        velVec = mvntVec(mvntThreshFrms:end);
        if max(~isnan(velVec)) == 0, continue, end
        [~,maxRef] = nanmax(velVec);
        if maxRef > numel(velVec)-winLeg, maxRef = numel(velVec)-winLeg;
        elseif maxRef < winLeg+1, maxRef = winLeg+1;
        end
        maxRef = maxRef(1);
        maxVel = velVec(maxRef);
        if maxVel > moveThresh, moveTest(iterV) = true;
        else nonMovers(iterV) = true;
        end
        if velVec(1) > stillThresh, earlyMovers(iterV) = true; end
        distCapVec = mvntArray(mvntThreshFrms+stimStart(iterV):end,1);
        distCapVec = distCapVec-distCapVec(1);
        distCapFrm = find(abs(distCapVec) > distCapThresh,1,'first');
%         if isempty(distCapFrm), distCapFrm = numel(distCapVec); end
        if isempty(distCapFrm)
            distCapFrm = NaN;
        end
        frmTotalOff = stimStart(iterV)+mvntThreshFrms-1;
        maxVelInfo(iterV,:) = [maxRef+frmTotalOff,...
            maxVel winLeg frmTotalOff distCapFrm+frmTotalOff];
    end
    
    earlyJumpers = fot < round(stimStart);
    earlyJumpers(isnan(earlyJumpers)) = false;
    if strcmp('Photoactivation',stimTypeString);
        earlyJumpers = fot < round(stimStart+frms_per_ms*10);
    end
    if strcmp(graphID,'noStim')
        data2keep = earlyJumpers;
    else
        jumpTest(earlyJumpers) = false;
        data2keep = jumpTest | nonJumpers;
    end
%     if graphOptions.keep_early_movers == 0
%         data2keep(~earlyMovers) = false;
%     end
    
    if sum(data2keep) == 0
        continue
    end
    graphTable.dataCount(graphID) = size(exptTable,1);
    exptTable = exptTable(data2keep,:);
    jumpTest = jumpTest(data2keep);
    moveTest = moveTest(data2keep);
    nonMovers = nonMovers(data2keep);
    earlyMovers = earlyMovers(data2keep);
    maxVelInfo = maxVelInfo(data2keep,:);
    shortLogical = shortLogical(data2keep);
    fowu_fowd_folp = fowu_fowd_folp(data2keep,:);
    fot = fot(data2keep);
    
    graphTable.videoList(graphID) = {exptTable.Properties.RowNames};
    graphTable.jumpTest(graphID) = {jumpTest};
    graphTable.shortTest(graphID) = {shortLogical};
    graphTable.fowu_fowd_folp(graphID) = {fowu_fowd_folp};
    graphTable.manCount(graphID) = sum(min(~isnan(fowu_fowd_folp),[],2));
    graphTable.jumpEnd(graphID) = {fot};
    graphTable.moveTest(graphID) = {moveTest & ~jumpTest};
    graphTable.nonMovers(graphID) = {nonMovers & ~jumpTest};
    graphTable.earlyMovers(graphID) = {earlyMovers};
    graphTable.maxVelInfo(graphID) = {maxVelInfo};
    graphTable.trkEnd(graphID) = {exptTable.trkEnd};
    graphTable.stimInfo(graphID) = {exptTable.stimInfo};
    graphTable.stimStart(graphID) = {exptTable.stimStart};
    graphTable.zeroFly_XYZmm_Tdeg_fac1000(graphID) = {exptTable.zeroFly_XYZmm_Tdeg_fac1000};
    graphTable.zeroFly_Trajectory(graphID) = {exptTable.zeroFly_Trajectory};
    graphTable.zeroFly_Jump(graphID) = {exptTable.zeroFly_Jump};
    graphTable.zeroFly_Departure(graphID) = {exptTable.zeroFly_Departure};
    graphTable.zeroFly_StimAtStimStart(graphID) = {exptTable.zeroFly_StimAtStimStart};
    graphTable.zeroFly_StimAtJump(graphID) = {exptTable.zeroFly_StimAtJump};
    graphTable.zeroFly_StimAtFrmOne(graphID) = {exptTable.zeroFly_StimAtFrmOne};
    graphTable.rawFly_XYZpix_Trad_frmOne(graphID) = {exptTable.rawFly_XYZpix_Trad_frmOne};
    graphTable.distFromCenter_pix(graphID) = {exptTable.distFromCenter_pix};
    graphTable.relMotion_FB_LR_Tdeg_fac100(graphID) = {exptTable.relMotion_FB_LR_Tdeg_fac100};
    graphTable.relPosition_FB_LR_Tdeg_fac100(graphID) = {exptTable.relPosition_FB_LR_Tdeg_fac100};
    graphTable.dist_speed_accel_fac100(graphID) = {exptTable.dist_speed_accel_fac100};
    graphTable.flyLength_mm(graphID) = {exptTable.flyLength_mm};
    graphTable.plotCount(graphID) = numel(jumpTest);
    graphTable.jumpCount(graphID) = sum(jumpTest);
    
    if strcmp('Visual_stimulation',stimTypeString) || strcmp('Combo',stimTypeString)
        if ~isempty(strfind(exptTable.stimInfo,'loom'))
            lvList = zeros(size(exptTable,1),1);
            for iterLV = 1:size(exptTable,1)
                stimInfo = exptTable.stimInfo{iterLV};
                stimInfo(1:strfind(stimInfo,'lv')+1) = [];
                stimInfo(strfind(stimInfo,'_'):end) = [];
                lvList(iterLV) = str2double(stimInfo);
            end
            graphTable.lvList(graphID) = {lvList};
        end
    end
end

[~,B] = xlsfinfo(excelPath);
sheetName = graphOptions.sheetName;
sheetValid = any(strcmp(B,sheetName));

if sheetValid
    write2Table = readtable(excelPath,'Sheet',sheetName,'ReadRowNames',true);
    savedNames = write2Table.Properties.RowNames;
    savedVars = write2Table.Properties.VariableNames;
    writeVars = {'vidCount','dataCount','plotCount','jumpCount','manCount'};
    writeVars = intersect(writeVars,savedVars);
    writeNames = intersect(graphTable.Properties.RowNames,savedNames);
    write2Table(writeNames,writeVars) = graphTable(writeNames,writeVars);
else
    writeVars = {'vidCount','dataCount','plotCount','jumpCount','manCount','recordRate',...
        'stimDur','genotype'};
    if max(strcmp(excelTable.Experiment_Type,'Visual_stimulation'))
        writeVars = cat(2,writeVars,{'elevation','azimuth','startSize','stopSize',...
            'lv','contrast','width','frequency','visStimType'});
    end
    write2Table = graphTable(:,writeVars);
    write2Table.plotID = write2Table.Properties.RowNames;
    write2Table.groupID = repmat({sheetName},size(write2Table,1),1);
    write2Table.order = (1:size(write2Table,1))';
end
try
    writetable(write2Table,excelPath,'Sheet',sheetName,'WriteRowNames',true)
catch
    disp('could not write to excel file - makeGraphingTable.mat')
end
graphTableRowNames = graphTable.Properties.RowNames;
bigNdx = cellfun(@(x) numel(x) > 63,graphTableRowNames);
graphTableRowNames(bigNdx) = cellfun(@(x) x(1:63),graphTableRowNames(bigNdx),'uniformoutput',false);
graphTable.Properties.RowNames = graphTableRowNames;

