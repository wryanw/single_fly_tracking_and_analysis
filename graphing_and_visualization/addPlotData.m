function [outputVar,excelTable] = addPlotData(graphTable)
%%
if ~exist('graphTable','var')
    graphTable = makeGraphingTable;
end
optionsPath = 'Z:\Data_pez3000_analyzed\WRW_graphing_variables\graphOptions.mat';
graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
excelPath = graphOptions.excelPath;
[~,excelName] = fileparts(excelPath);
sheetName2plot = graphOptions.sheetName;
excelTable = readtable(excelPath,'Sheet',sheetName2plot,'ReadRowNames',true);
dataIDlist = excelTable.Properties.RowNames;
if max(strcmp(excelTable.Properties.VariableNames,'keep'))
    dataIDlist = dataIDlist(logical(excelTable.keep));
end
dataIDlist = intersect(dataIDlist,graphTable.Properties.RowNames);
plotTable = graphTable(dataIDlist,:);
graphOptions.dataIDlist = dataIDlist;
dataCt = size(dataIDlist,1);
if dataCt == 0
    error('zero names from excel file match row names in graphTable')
end
graphName = graphOptions.graphName;
subset2plot = graphOptions.subset2plot;
jump_thresh = graphOptions.jump_thresh;
keep_late_jumpers = graphOptions.keep_late_jumpers;
if isfield(graphOptions,'padOption')
    padOption = graphOptions.padOption;
else
    padOption = 1;
end
deg2rad = @(x) x*(pi/180);
rad2deg = @(x) x*(180/pi);
% medDur = median(plotTable.stimDur);
plotTable.jumpCount = zeros(size(plotTable,1),1);
plotTable.shortUseTotal = zeros(size(plotTable,1),1);
plotTable.returnData = cell(size(plotTable,1),1);
graphTableRowNames = plotTable.Properties.RowNames;
bigNdx = cellfun(@(x) numel(x) > 63,graphTableRowNames);
graphTableRowNames(bigNdx) = cellfun(@(x) x(1:63),graphTableRowNames(bigNdx),'uniformoutput',false);
plotTable.Properties.RowNames = graphTableRowNames;
% [~,excelName] = fileparts(excelPath);
% medFlyLength = median(cell2mat(plotTable.flyLength_mm(dataIDlist)));
subRefList = {'videoList'
    'jumpTest'
    'shortTest'
    'moveTest'
    'nonMovers'
    'earlyMovers'
    'maxVelInfo'
    'fowu_fowd_folp'
    'jumpEnd'
    'trkEnd'
    'stimInfo'
    'stimStart'
    'zeroFly_XYZmm_Tdeg_fac1000'
    'zeroFly_Trajectory'
    'zeroFly_Jump'
    'zeroFly_Departure'
    'zeroFly_StimAtStimStart'
    'zeroFly_StimAtJump'
    'zeroFly_StimAtFrmOne'
    'rawFly_XYZpix_Trad_frmOne'
    'dist_speed_accel_fac100'
    'distFromCenter_pix'
    'relMotion_FB_LR_Tdeg_fac100'
    'relPosition_FB_LR_Tdeg_fac100'
    'lvList'
    'flyLength_mm'};
for iterD = 1:dataCt
    dataID = dataIDlist{iterD};
    if plotTable.plotCount(dataID) == 0
        continue
    end
    stimTypeString = plotTable.stimType{dataID};
    
    if strcmp('Visual_stimulation',stimTypeString)
        exptAzi = plotTable.azimuth(dataID);
        exptEle = plotTable.elevation(dataID);
    else
        exptAzi = 0;
        exptEle = 90;
    end
    
    frms_per_ms = round(plotTable.recordRate(dataID)/1000);
    frm2ms = @(x) x/frms_per_ms;
    ms2frm = @(x) round(x*frms_per_ms);
    jumpTest = plotTable.jumpTest{dataID};
    nonJumpers = ~jumpTest;
    moveTest = plotTable.moveTest{dataID};
    nonMovers = plotTable.nonMovers{dataID};
    earlyMovers = plotTable.earlyMovers{dataID};
    shortLogical = plotTable.shortTest{dataID};
    plotTable.jumpCount(dataID) = sum(jumpTest);
    plotTable.shortUseTotal(dataID) = sum(shortLogical(jumpTest));
    %%%%% scatter plot stuff
    if strcmp(subset2plot,'jumping')
        data2keep = jumpTest;
    elseif strcmp(subset2plot,'nonjumping')
        data2keep = nonJumpers;
    elseif strcmp(subset2plot,'moving')
        data2keep = moveTest;
    elseif strcmp(subset2plot,'nonmoving')
        data2keep = nonMovers;
    else
        data2keep = jumpTest | nonJumpers;
    end
    fot = plotTable.jumpEnd{dataID};
    fowu_fowd_folp = plotTable.fowu_fowd_folp{dataID};
    fowu = fowu_fowd_folp(:,1);
    fowd = fowu_fowd_folp(:,2);
    folp = fowu_fowd_folp(:,3);
    maxVelInfo = plotTable.maxVelInfo{dataID};
    
    stimInit = plotTable.zeroFly_StimAtStimStart{dataID};
    stimInit = centerDegs(stimInit);
    stimStart = plotTable.stimStart{dataID};
    
%     stimDur = plotTable.stimDur(dataID);
    stimDur = excelTable.stimDur(dataID);
    
    lateJumpers = false(size(jumpTest));
    if strcmp('Photoactivation',stimTypeString);
        lateJumpers = fot > round(stimStart+jump_thresh*frms_per_ms);
        lateJumpers(isnan(lateJumpers)) = false;
    elseif strcmp(excelName,'DL_looming')
        jump_thresh = 455;
        lateJumpers = fot > round(stimStart+jump_thresh*frms_per_ms);
        lateJumpers(isnan(lateJumpers)) = false;
    elseif ~keep_late_jumpers
        lateJumpers = fot > round(stimStart+stimDur);
        lateJumpers(isnan(lateJumpers)) = false;
    end
    data2keep(lateJumpers) = false;
%     if ~isinf(jump_thresh)
%         stimDur_ms = jump_thresh;
%         stimDur = zeros(size(stimDur))+ms2frm(stimDur_ms);
%     end
    if graphOptions.keep_early_movers == 0
        data2keep(earlyMovers) = false;
    end
    if ~isempty(strfind(graphName,'jaabaReaching'))
        exptIDlist = plotTable.exptIDs{dataID};
        returnData = cell(numel(exptIDlist),1);
        trkEnd = round(plotTable.trkEnd{dataID}/10);
        for iterJ = 1:numel(exptIDlist)
            exptID = exptIDlist{iterJ}(4:end);
            jaabaTablePath = ['Z:\Data_pez3000_analyzed\' exptID '\' exptID '_JAABAannotations.mat'];
            if ~exist(jaabaTablePath,'file')
                disp([exptID ' - no jaaba file'])
                continue
            end
            jaabaTableLoading = load(jaabaTablePath);
            try
                jaabaTableRaw = jaabaTableLoading.saveObj;
            catch
                jaabaTableRaw = jaabaTableLoading.jaabaTable;
            end
            jaabaVids = jaabaTableRaw.Properties.RowNames;
            vids2keep = plotTable.videoList{dataID};
            jaabaData = cell(numel(vids2keep),1);
            for iterV = 1:numel(vids2keep)
                videoID = vids2keep{iterV};
                if ~max(strcmp(jaabaVids,videoID)), continue, end
                jaabaVec = jaabaTableRaw.Rearing(videoID).postprocessed{1};
                jaabaVec(trkEnd(iterV):end) = NaN;
                jaabaVec = jaabaVec(round(stimStart(iterV)/10):end);
                jaabaData{iterV} = jaabaVec;
            end
            returnData = cat(1,returnData,jaabaData);
        end
        maxL = max(cellfun(@(x) numel(x),returnData));
        for iterR = 1:numel(returnData)
            if isempty(returnData{iterR}), continue, end
            nanPad = NaN(1,maxL-numel(returnData{iterR})+1);
            returnData{iterR} = cat(2,returnData{iterR},nanPad);
        end
        returnData = cat(1,returnData{:});
        returnData(:,sum(~isnan(returnData)) < 10) = NaN;
        returnData = nanmean(returnData);
        returnData = repmat(returnData,10,1);
        returnData = returnData(:)';
    elseif ~isempty(strfind(graphName,'initialOrientation'))
        returnData = plotTable.rawFly_XYZpix_Trad_frmOne{dataID}(data2keep,4);
    elseif ~isempty(strfind(graphName,'flyLength'))
        returnData = plotTable.flyLength_mm{dataID}(data2keep);
    elseif strcmp(graphName(1:6),'escape')
        nameParts = strsplit(graphName,'_');
        if isempty(strfind(graphName,'Meatball')) %%%% negative test
            escData = plotTable.(['zeroFly_' regexprep(nameParts{1},'escape','')]){dataID};
            zeroFlyVec = escData(:,1);
            eleFlyVec = escData(:,3);
            rhoVec = escData(:,2);
            eleRhoVec = escData(:,4);
        elseif ~isempty(strfind(graphName,'Meatball'))
            posCell = plotTable.zeroFly_XYZmm_Tdeg_fac1000{dataID};
            locStimStart = stimStart;
            zeroFlyVec = zeros(numel(posCell),1);
            rhoVec = zeros(numel(posCell),1);
            eleRhoVec = zeros(numel(posCell),1);
            eleFlyVec = zeros(numel(posCell),1);
            parfor iterR = 1:numel(posCell)
                data2keep(iterR) = false;
                mvntArray = double(posCell{iterR}(:,1:3))/1000;
                if strcmp('Photoactivation',stimTypeString)
                    dataBegin = locStimStart(iterR);
                else
                    dataBegin = maxVelInfo(iterR,4);
                end
                dataBegin = round(dataBegin);
                if isnan(dataBegin), continue, end
                if dataBegin < 1, continue, end
                if dataBegin > size(mvntArray,1), continue, end
                distCapVal = maxVelInfo(iterR,5);
                if distCapVal <= dataBegin, continue, end
                if ~isnan(distCapVal)
                    mvntArray(distCapVal:end,:) = [];
                end
                dataPts = mvntArray(end,:)-mvntArray(dataBegin,:);
                [theta,ele,elerho] = cart2sph(dataPts(1),dataPts(2),dataPts(3));
                zeroFlyVec(iterR) = rad2deg(theta);
                eleFlyVec(iterR) = rad2deg(ele);
                [~,rho] = cart2pol(dataPts(1),dataPts(2));
                rhoVec(iterR) = rho;
                eleRhoVec(iterR) = elerho;
                data2keep(iterR) = true;
            end
        end
        if strcmp(nameParts{2},'elevation')
            returnData = [zeroFlyVec(data2keep) rhoVec(data2keep) eleFlyVec(data2keep) eleRhoVec(data2keep)];
        else
            if strcmp(nameParts{2},'zeroStimAzi') && ~strcmp('Photoactivation',stimTypeString)
                dirResult = zeroFlyVec+stimInit-exptAzi;
            elseif strcmp(nameParts{2},'3D')
                [x1,y1,z1] = sph2cart(deg2rad(exptAzi),deg2rad(exptEle),1);
                dirResult = NaN(size(zeroFlyVec,1),1);
                
                %%%%% For modeling mechanical bias only !!! %%%%%
                if ~isempty(strfind(graphName,'_NoElev')) || ~isempty(strfind(graphName,'_AllBiophys'))
                    eleFlyVec(data2keep) = 50;
                end
                if ~isempty(strfind(graphName,'_NoAzim')) || ~isempty(strfind(graphName,'_AllBiophys'))
                    zeroFlyVec(data2keep) = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for iterT = 1:size(zeroFlyVec,1)
                    [x2,y2,z2] = sph2cart(deg2rad(zeroFlyVec(iterT)),...
                        deg2rad(eleFlyVec(iterT)),1);
                    a = [x1,y1,z1];
                    b = [x2,y2,z2];
                    dirResult(iterT) = atan2(norm(cross(a,b)),dot(a,b));
                end
                dirResult = rad2deg(dirResult);
            else
                dirResult = zeroFlyVec;
            end
            data2keep(isnan(dirResult)) = false;
            dirResult = centerDegs(dirResult);
            returnData = [dirResult(data2keep) rhoVec(data2keep) eleFlyVec(data2keep) eleRhoVec(data2keep)];
            if strcmp(subset2plot,'jumping')
                plotTable.jumpCount(dataID) = size(returnData,1);
            else
                plotTable.plotCount(dataID) = size(returnData,1);
            end
        end
    elseif max(cellfun(@(x) ~isempty(strfind(graphName,x)),{'relMotion','absMotion','spaghettiPlot'}))
        lengthArray = plotTable.flyLength_mm{dataID}(data2keep);
        visStimType = excelTable.visStimType{dataID};
        if ~isempty(strfind(graphName,'spaghettiPlot'))
            posCell = plotTable.zeroFly_XYZmm_Tdeg_fac1000{dataID}(data2keep);
        elseif ~isempty(strfind(graphName,'relMotion')) || ~isempty(strfind(graphName,'absMotion_turning'))
%             if strcmp(visStimType,'grating')
%                 posCell = plotTable.relMotion_FB_LR_Tdeg_fac100{dataID}(data2keep);
%                 turnCell = plotTable.zeroFly_XYZmm_Tdeg_fac1000{dataID}(data2keep);
%             else
            if ~isempty(strfind(graphName,'Dist')) || ~isempty(strfind(graphName,'absMotion_turning_Vel'))
                posCell = plotTable.relPosition_FB_LR_Tdeg_fac100{dataID}(data2keep);
            else
                posCell = plotTable.relMotion_FB_LR_Tdeg_fac100{dataID}(data2keep);
            end
        else
            posCell = plotTable.dist_speed_accel_fac100{dataID}(data2keep);
        end
        locFot = fot(data2keep);
        fowu_fowd_folp = fowu_fowd_folp(data2keep,:);
        locStimStart = stimStart(data2keep);
        dataCell = cell(numel(posCell),1);
        subset2keep = data2keep(data2keep);
        useManual = graphOptions.useManual;
        initStimSize = plotTable.startSize(dataID);
        %             ellovervee = graphTable.lvList{dataID}(1);
        ellovervee = plotTable.lv(dataID);
        finalStimSize = plotTable.stopSize(dataID);
        maxVelInfo = maxVelInfo(data2keep,:);
        for iterR = 1:numel(posCell)
            subset2keep(iterR) = false;
            manTest = ~max(isnan(fowu_fowd_folp(iterR,:))) && useManual ~= 2;
            if ~isempty(strfind(graphName,'spaghettiPlot'))
                mvntArray = double(posCell{iterR}(:,1:2))/1000;
            else
                mvntArray = double(posCell{iterR})/100;
%                 turnArray = double(turnCell{iterR});
            end
            dataBegin = locStimStart(iterR);
            dataEnd = locStimStart(iterR)+stimDur;
            if ~isempty(strfind(graphName,'timeBounds'))
                strA = graphName(strfind(graphName,'timeBounds')+10:end);
                dataBegin = ms2frm(str2double(strA(1:strfind(strA,'to')-1)));
                strB = strA(strfind(strA,'to')+2:end);
                if isempty(strfind(strB,'_'))
                    dataEnd = ms2frm(str2double(strB(1:end)));
                else
                    dataEnd = ms2frm(str2double(strB(1:strfind(strB,'_')-1)));
                end
            elseif ~isempty(strfind(graphName,'MidLoom')) && ~strcmp('Photoactivation',stimTypeString)
                dataBegin = locStimStart(iterR)+ms2frm(0);
                dataEnd = dataBegin+ms2frm(250);
            elseif ~isempty(strfind(graphName,'FullDur')) && ~strcmp('Photoactivation',stimTypeString)
                dataEnd = size(mvntArray,1);
            elseif ~isempty(strfind(graphName,'LateLoom')) && ~strcmp('Photoactivation',stimTypeString)
                dataBegin = maxVelInfo(iterR,4);
            elseif ~isempty(strfind(graphName,'Deg2deg')) % currently not used
                stimTimeStep = (1/360)*1000;%milliseconds per frame channel at 120 Hz
                minTheta = deg2rad(initStimSize);
                maxTheta = deg2rad(finalStimSize);
                stimStartTime = ellovervee/tan(minTheta/2);
                stimEndTime = ellovervee/tan(maxTheta/2);
                stimTimeVector = fliplr(stimEndTime:stimTimeStep:stimStartTime);
                stimThetaVector = 2*atan(ellovervee./stimTimeVector)./(pi/180);
                stimTimeVector = fliplr(stimTimeVector);
                tagNdx = strfind(graphName,'Deg2deg');
                degStart = str2double(graphName(tagNdx+8:tagNdx+10));
                degStop = str2double(graphName(tagNdx+12:tagNdx+14));
                dataBegin = ms2frm(round(stimTimeVector(find(stimThetaVector > degStart,1,'first'))))+stimStart;
                dataEnd = ms2frm(round(stimTimeVector(find(stimThetaVector < degStop,1,'last'))))+stimStart;
            elseif ~isempty(strfind(graphName,'FullJumpDur'))
                dataBegin = locFot(iterR)-round(frms_per_ms*9);
                dataEnd = locFot(iterR)+round(frms_per_ms*1);
            elseif ~isempty(strfind(graphName,'Fot'))
                meanSpan = frms_per_ms/2;
                dataBegin = locFot(iterR)-meanSpan/2;
                dataEnd = dataBegin+frms_per_ms;
            elseif ~isempty(strfind(graphName,'MidLegExt'))
                if manTest
                    frameList = linspace(fowu_fowd_folp(iterR,3),locFot(iterR),4);
                    dataBegin = frameList(2);
                    dataEnd = frameList(3);
                else
                    dataBegin = locFot(iterR)-frms_per_ms*4.5;
                    dataEnd = dataBegin+frms_per_ms;
                end
            elseif ~isempty(strfind(graphName,'FullLegExt'))
                if manTest
                    dataBegin = fowu_fowd_folp(iterR,3);
                    dataEnd = locFot(iterR);
                else
                    dataBegin = locFot(iterR)-frms_per_ms*4.5;
                    dataEnd = locFot(iterR);
                    if dataEnd > size(mvntArray,1)
                        dataEnd = size(mvntArray,1);
                        dataBegin = dataEnd-frms_per_ms*4.5;
                    end
                end
            elseif ~isempty(strfind(graphName,'MaxVelWin'))
                dataBegin = maxVelInfo(iterR,1)-maxVelInfo(iterR,3);
                dataEnd = maxVelInfo(iterR,1)+maxVelInfo(iterR,3);
            end
            if ~isempty(strfind(graphName,'JumpCap')) && ~isnan(locFot(iterR))
                if manTest
                    dataEnd = fowu_fowd_folp(iterR,3);
                else
                    dataEnd = locFot(iterR)-frms_per_ms*10;
                end
            end
            dataBegin = round(dataBegin);
            dataEnd = round(dataEnd);
            dataLength = dataEnd-dataBegin+1;
            
            if isnan(dataBegin), continue, end
            if dataBegin < 1, continue, end
            if dataLength < 1, continue, end
            if size(mvntArray,1) < 10, continue, end
            if dataBegin > size(mvntArray,1), continue, end
            if ~isempty(strfind(graphName,'MvntCap'))
                distCapVal = maxVelInfo(iterR,5);
                if distCapVal <= dataBegin, continue, end
                if ~isnan(distCapVal)
                    mvntArray(distCapVal:end,:) = [];
                end
            end
            if ~isempty(strfind(graphName,'spaghettiPlot'))
                mvntArray(isnan(mvntArray(:,1)),:) = [];
                if dataEnd > size(mvntArray,1) || isnan(dataEnd)
                    dataEnd = size(mvntArray,1);
                end
                dataArray = mvntArray(dataBegin:dataEnd,:);
                dataArray = dataArray-repmat(dataArray(1,:),size(dataArray,1),1);
                dataCell(iterR) = {cat(1,dataArray,[NaN NaN])};
            elseif ~isempty(strfind(graphName,'relMotion')) || ~isempty(strfind(graphName,'absMotion_turning'))
                if ~isempty(strfind(graphName,'walking'))
                    dataVec = mvntArray(dataBegin:end,1);
                elseif ~isempty(strfind(graphName,'movement'))
                    dataVec = sqrt(sum(mvntArray(dataBegin:end,1:2).^2,2));
                elseif ~isempty(strfind(graphName,'sidepass'))
                    dataVec = mvntArray(dataBegin:end,2);
                    dataVec = dataVec*(-1);%this makes negative mean 'away from stimulus'
                elseif ~isempty(strfind(graphName,'turning'))
                    dataVec = mvntArray(dataBegin:end,3);
                    if ~isempty(strfind(graphName,'Vel_'))
%                         dataVec = dataVec/1000;
                    end
                end
                if strcmp('Photoactivation',stimTypeString) && (~isempty(strfind(graphName,'sidepass')) || ~isempty(strfind(graphName,'turning')))
                    dataVec = abs(dataVec);
                end
                if ~isempty(strfind(graphName,'absMotion_turning'))
                    dataVec = abs(dataVec);
                end
                if dataLength > numel(dataVec)
                    dataLength = numel(dataVec);
                end
                dataVec = dataVec(1:dataLength);
                if ~isempty(strfind(graphName,'vector'))
                    if ~isempty(strfind(graphName,'Dist'))
                        dataVec = dataVec-dataVec(1);
                    end
                    dataCell(iterR) = {dataVec(:)'};
                else
                    dataVec(isnan(dataVec)) = 0;
                    if ~isempty(strfind(graphName,'absMotion_turning_Vel'))
                        posVal = abs(dataVec(end)-dataVec(1)/numel(dataVec));
                    elseif ~isempty(strfind(graphName,'Vel_'))
                        posVal = mean(dataVec);
                    else
                        posVal = (dataVec(end)-dataVec(1));
                    end
                    if (~isempty(strfind(graphName,'sidepass')) || ~isempty(strfind(graphName,'turning'))) && ~strcmp(visStimType,'grating')
                        posVal = abs(posVal);
                    end
                    dataCell(iterR) = {posVal};
                end
            else
                if ~isempty(strfind(graphName,'Accel')) || ~isempty(strfind(graphName,'Force'))
                    dataVec = mvntArray(dataBegin:end,3);
                elseif ~isempty(strfind(graphName,'Vel_'))
                    dataVec = mvntArray(dataBegin:end,2);
                else
                    dataVec = mvntArray(dataBegin:end,1);
                    dataVec = (dataVec-dataVec(1));
                end
                if ~isempty(strfind(graphName,'Force'))
                    dataVec = dataVec*lengthArray(iterR)^2/5.78*0.89;% micro Newtons
                end
                if ~isempty(strfind(graphName,'vector'))
                    if dataLength <= numel(dataVec)
                        dataVec = dataVec(1:dataLength);
                    end
                else
                    if dataLength > numel(dataVec)
                        dataLength = numel(dataVec);
                    end
                    if ~isempty(strfind(graphName,'Dist'))
                        dataVec = dataVec(dataLength);%total distance mm
                    else
                        dataVec = mean(dataVec(1:dataLength));
                    end
                end
                if ~isempty(dataVec)
                    if ~isempty(strfind(graphName,'FullDur'))
                        dataL = numel(dataVec);
                        if dataL >= stimDur
                            dataVec = dataVec(1:stimDur);
                        end
                    end
                end
                dataCell(iterR) = {dataVec(:)'};
            end
            subset2keep(iterR) = true;
        end
        if ~isempty(strfind(graphName,'vector'))
            dataLength = max(cellfun(@(x) numel(x),dataCell));
            for iterL = 1:numel(dataCell)
                dataVec = dataCell{iterL};
                if dataLength > numel(dataVec)
                    if isempty(dataVec), continue, end
                    if padOption == 1
                        dataVec = cat(1,dataVec(:),NaN(dataLength-numel(dataVec),1));
                    else
                        dataVec = cat(1,NaN(dataLength-numel(dataVec),1),dataVec(:));
                    end
                end
                dataCell(iterL) = {dataVec(:)'};
            end
        end
        returnData = cat(1,dataCell{:});
        if ~isempty(strfind(graphName,'spaghettiPlot'))
            plotTable.plotCount(dataID) = sum(isnan(returnData(:,1)));
        elseif strcmp(subset2plot,'jumping')
            plotTable.jumpCount(dataID) = size(returnData,1);
        else
            plotTable.plotCount(dataID) = size(returnData,1);
        end
    elseif ~isempty(strfind(graphName,'departureScatter'))
        if strcmp('Photoactivation',stimTypeString); continue, end
        zeroFlyVec = plotTable.zeroFly_Jump{dataID}(:,1);
        stimJump = plotTable.zeroFly_StimAtJump{dataID};
        dirResult = zeroFlyVec-(stimJump);
        while sum(dirResult > 0) > 0,dirResult(dirResult > 0) = dirResult(dirResult > 0)-360;end
        while sum(dirResult < -360) > 0,dirResult(dirResult < -360) = dirResult(dirResult < -360)+360;end
        dirResult = dirResult+(stimJump);
        data2keep(isnan(dirResult)) = false;
        returnData = [dirResult(data2keep) stimJump(data2keep)];
        if strcmp(subset2plot,'jumping')
            plotTable.jumpCount(dataID) = size(returnData,1);
        else
            plotTable.plotCount(dataID) = size(returnData,1);
        end
    elseif ~isempty(strfind(graphName,'frameOf')) || ~isempty(strfind(graphName,'durationOf'))
        normalFrm = [fowu folp fowd fot maxVelInfo(:,1) maxVelInfo(:,5)]-repmat(stimStart,1,6);
        if strcmp(subset2plot,'jumping')
%             normalFrm(:,5) = normalFrm(:,4)-normalFrm(:,5);
        end
%         normalFrm = normalFrm-ms2frm(8.5); %%%%%%%%%seems to be a 8.5 ms delay
        normalFrm = frm2ms(normalFrm);
        if ~isempty(strfind(graphName,'time')) || ~isempty(strfind(graphName,'durationOf'))
            normTo180fot = normalFrm;
            if strcmp('Visual_stimulation',stimTypeString) || strcmp('Combo',stimTypeString)
                if strcmp(plotTable.visStimType{dataID},'loom')
                    initStimSize = plotTable.startSize(dataID);
                    if strcmp('Combo',stimTypeString)
                        ellovervee = zeros(size(normTo180fot,1),1)+plotTable.lv(dataID);
                    else
                        ellovervee = plotTable.lvList{dataID};
                    end
                    minTheta = deg2rad(initStimSize);
                    tocTime = ellovervee/tan(minTheta/2);
                    normTo180fot = normTo180fot-repmat(tocTime,1,6);
                end
            else
%                 if ~exist('actKey','var')
%                     actKey = readtable(excelPath,'ReadRowNames',false,'Sheet','activation_key');
%                     actKey = table2cell(actKey);
%                     actKey = actKey(:,1:2);
%                     varDir = [filesep filesep 'DM11' filesep 'cardlab' filesep,...
%                         'pez3000_variables' filesep 'photoactivation_stimuli'];
%                 end
%                 if  graphOptions.sliceOption == 1
%                     var2use = strcmp(actKey(:,2),plotTable.stimInfo{dataID}{1});
%                     varName = actKey{var2use,1};
%                 else
%                     varName = plotTable.stimInfo{dataID}{1};
%                 end
%                 if ~isempty(strfind(varName,'ramp'))
%                     varPath = fullfile(varDir,[varName '.mat']);
%                     photoVar = load(varPath);
%                     varDur = ms2frm(photoVar.var_tot_dur*1.028);
%                     prctPwr = zeros(1,varDur);
%                     rampVals = linspace(photoVar.var_ramp_init,photoVar.var_intensity,ms2frm(photoVar.var_ramp_width*1.028));
%                     prctPwr(1:numel(rampVals)) = rampVals;
%                     stimStartOffset = find(prctPwr >= 1,1,'first');
%                     normTo180fot = normTo180fot-stimStartOffset;
%                 end
            end
        else
            if strcmp('Photoactivation',stimTypeString);
                error('incompatible stimulus')
            end
            initStimSize = plotTable.startSize(dataID);
            ellovervee = plotTable.lvList{dataID};
            minTheta = deg2rad(initStimSize);
            stimStartTime = ellovervee/tan(minTheta/2);
            fotTime = repmat(stimStartTime,1,6)-normalFrm;
            fotTheta = 2.*atan(repmat(ellovervee,1,6)./fotTime);
            normTo180fot = rad2deg(fotTheta);
            normTo180fot(fotTime < 0) = normTo180fot(fotTime < 0)+360;
        end
        returnData = normTo180fot(data2keep,:);
    elseif ~isempty(strfind(graphName,'jumpingJoules')) || ~isempty(strfind(graphName,'jumpingWatts'))
        if ~strcmp('Photoactivation',stimTypeString);
            error('incompatible stimulus')
        end
        normalFrm = [fot maxVelInfo(:,1) maxVelInfo(:,5)]-repmat(stimStart,1,3);
        if ~exist('actKey','var')
            actKey = readtable(excelPath,'ReadRowNames',false,'Sheet','activation_key');
            actKey = table2cell(actKey);
            actKey = actKey(:,1:2);
            varDir = [filesep filesep 'DM11' filesep 'cardlab' filesep,...
                'pez3000_variables' filesep 'photoactivation_stimuli'];
        end
        if ~exist('fitFun','var')
            dataDir = 'Z:\pez3000_variables';
            fitFunLoad = load(fullfile(dataDir,'photoactivationTransferFunction.mat'),'fitFun');
            fitFun = fitFunLoad.fitFun;
            funStr = func2str(fitFun);
            funStr = cat(2,funStr(1:4),'4*(',funStr(5:end),')');
            fitFun = str2func(funStr);
        end
        prct2mW = fitFun;% generated from avg measurements of LEDs
        if  graphOptions.sliceOption == 1
            var2use = strcmp(actKey(:,2),plotTable.stimInfo{dataID}{1});
            varName = actKey{var2use,1};
        else
            varName = plotTable.stimInfo{dataID}{1};
        end
        varPath = fullfile(varDir,[varName '.mat']);
        photoVar = load(varPath);
        if ~isempty(strfind(varName,'ramp'))
            varDur = ms2frm(photoVar.var_tot_dur*1.028);
            prctPwr = zeros(1,varDur);
            rampVals = linspace(photoVar.var_ramp_init,photoVar.var_intensity,ms2frm(photoVar.var_ramp_width*1.028));
            prctPwr(1:numel(rampVals)) = rampVals;
            prctPwr(prctPwr < 3) = 0;
            Woutput = prct2mW(prctPwr);
            Woutput(Woutput < 0) = 0;
            Woutput = Woutput/1000;% converts from uW to mW
            power = Woutput;
        elseif ~isempty(strfind(varName,'pulse'))
            varDur = ms2frm(photoVar.var_tot_dur*1.028);
            pulsePos = photoVar.zPulse;
            prctPwr = photoVar.aPulse;
            Woutput = prct2mW(prctPwr);
            Woutput(Woutput < 0) = 0;
            Woutput = Woutput/1000;% converts from uW to mW
            power = zeros(1,varDur)+Woutput(end);
            pulsePos(1) = 1;
            for i = 1:numel(pulsePos)-1
                startPt = ms2frm(pulsePos(i)*1.028);
                endPt = ms2frm(pulsePos(i+1)*1.028);
                power(startPt:endPt) = Woutput(i);
            end
        else
            error('unknown photoactivation stimulus')
        end
        if ~isempty(strfind(graphName,'jumpingJoules'))
            powerCum = cumsum(power);%integrate
            powerCum = frm2ms(powerCum);%necessary because of integration
            powerCum = powerCum/1000;%convert to sec from ms since Joule*sec is Watt
        else
            powerCum = power;
        end
        for iterRa = 1:size(normalFrm,1)
            for iterRb = 1:size(normalFrm,2)
                if isnan(normalFrm(iterRa,iterRb))
                    continue
                end
                if normalFrm(iterRa,iterRb) > numel(powerCum)
                    normalFrm(iterRa,iterRb) = powerCum(end);
                elseif normalFrm(iterRa,iterRb) < 1
                    normalFrm(iterRa,iterRb) = powerCum(1);
                else
                    normalFrm(iterRa,iterRb) = powerCum(normalFrm(iterRa,iterRb));
                end
            end
        end
        
        returnData = normalFrm(data2keep,:);
    elseif ~isempty(strfind(graphName,'maxMotion'))
        returnData = maxVelInfo(data2keep,2);
    elseif ~isempty(strfind(graphName,'deltaHeading'))
        deltaHeadingNet = abs(plotTable.zeroFly_StimAtJump{dataID}-plotTable.zeroFly_StimAtStimStart{dataID});
        returnData = deltaHeadingNet(data2keep);
    elseif ~isempty(strfind(graphName,'pitch'))
        returnData = plotTable.pitch75frm(data2keep,1);
    elseif ~isempty(strfind(graphName,'fromCenter'))
        returnData = plotTable.distFromCenter_pix(data2keep);
    elseif strcmp(graphName,'wing2fot')
        %             x_start = 0;    x_end = 7;  x_step = .075;
        histVals = hist(histVals,(x_start:x_step:x_end));
        returnData = fot(data2keep)-fowu(data2keep);
        %             returnData = histVals+rand(size(histVals));
    elseif ~isempty(strfind(graphName,'wing2fot_histogram'))
        histVals = (fot(data2keep)-fowu(data2keep));
%         histVals = (fowd(data2keep)-fowu(data2keep));
        histVals = real(log(histVals./6));
        if ~isempty(strfind(graphName,'_VSazi'))
            returnData = [histVals stimInit(data2keep)];
        else
            returnData = histVals;
        end
    elseif strcmp(graphName,'groundTruthing_all')
        returnData = [plotTable.centerDist(data2keep) fot(data2keep),...
            plotTable.maxDist(data2keep) plotTable.trkEnd(data2keep)];
    else % has info for montage
        returnArray = [fowu(:) fot stimStart repmat(stimDur,size(fot,1),1) plotTable.trkEnd{dataID}(:) fowd(:) folp(:)];
        returnData = returnArray(data2keep,:);
    end
    plotTable.returnData(dataID) = {returnData};
    for iterS = 1:numel(subRefList)
        if strcmp('Photoactivation',stimTypeString) && strcmp(subRefList{iterS},'lvList')
            continue
        end
        subArray = plotTable.(subRefList{iterS}){dataID}(data2keep,:);
        if exist('subset2keep','var')
            subArray = subArray(subset2keep,:);
        end
        plotTable.(subRefList{iterS}){dataID} = subArray;
    end
end

plotTable.groupID(dataIDlist) = excelTable.groupID(dataIDlist);
% plotTable.plotID(dataIDlist) = matlab.lang.makeUniqueStrings(excelTable.plotID(dataIDlist));
plotTable.plotID(dataIDlist) = excelTable.plotID(dataIDlist);
plotTable.order(dataIDlist) = excelTable.order(dataIDlist);
if max(strcmp(excelTable.Properties.VariableNames,'repeat'))
    plotTable.repeat(dataIDlist) = logical(excelTable.repeat(dataIDlist) == 1);
end

plotIDlist = unique(plotTable.plotID);
newOrder = zeros(numel(plotIDlist),1);
for iterU = 1:numel(plotIDlist)
    orderNdx = find(strcmp(plotTable.plotID,plotIDlist{iterU}),1,'first');
    newOrder(iterU) = plotTable.order(orderNdx);
end
[~,sortNdcs] = sort(newOrder);
newLabels = plotIDlist(sortNdcs);
newOrder = (1:numel(newLabels));

groupIDlist = unique(plotTable.groupID);
for iterG = 1:numel(groupIDlist)
    grpID = groupIDlist{iterG};
    grpBool = strcmp(grpID,plotTable.groupID);
    grpTable = plotTable(grpBool,:);
    grpIDlist = grpTable.plotID;
    [Lia,Lib] = ismember(newLabels,grpIDlist);
    Lib(Lib == 0) = [];
    grpTable.order(Lib) = newOrder(Lia);
    plotTable(grpBool,:) = grpTable;
end
if max(strcmp(excelTable.Properties.VariableNames,'repeat'))
    if sum(plotTable.repeat == 2) > 1
        error('too many plot repeats; can be no more than one')
    end
    if sum(excelTable.repeat == 2) == 1
        name2repeat = excelTable.Properties.RowNames{excelTable.repeat == 2};
        table2repeat = plotTable(name2repeat,:);
        plotTable(name2repeat,:) = [];
        names2add = matlab.lang.makeUniqueStrings(repmat({name2repeat},numel(plotIDlist),1));
        tables2add = cell(numel(names2add),1);
        for iterA = 1:numel(names2add)
            table2repeat.Properties.RowNames = names2add(iterA);
            tables2add{iterA} = table2repeat;
        end
        tables2add = cat(1,tables2add{:});
        tables2add.plotID = plotIDlist(sortNdcs);
        tables2add.order = newOrder(:);
        plotTable = cat(1,plotTable,tables2add);
    end
end
graphOptions.plotIDlist = plotIDlist(sortNdcs);
graphOptions.groupIDlist = groupIDlist;
save(optionsPath,'graphOptions')

if isempty(mfilename) || nargin == 0
    outputVar = [];
else
    outputVar = plotTable;
end