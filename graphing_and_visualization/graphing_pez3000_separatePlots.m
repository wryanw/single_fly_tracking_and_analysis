function pdfName = graphing_pez3000_separatePlots
%%
graphTable = makeGraphingTable;
%%
clearvars -except graphTable
rng('default');
rng(19);

segmentChoice = 8; % 1-full 10 ms...2-fot...3-mid leg extend...4-full leg
% extend...5-entire vector...6-half loom...7-jumpCap...8-last half
% loom...9-distCap

% use with escape plot
dirChoice = [1 4 1 1]; % [1-Jump 2-Departure 3-Trajectory 4-Meatball    1-zeroStimAzi 2-zeroFlyAzi 3-3D 4-elevation
%    1-Angle 2-Dist 3-Binomial    1-Variable Rho 2-One Rho]

plotChoice = 3;
groupChoice = 1;
saveFig = 0;
singleAx = 0;
rhoThresh = 0.0;%0.05
plotCircHist = 0;
plotCircArrow = 0;
circConfStyle = 1;
showDotsOrLines = 1;
gaussFitCt = 3;
barwidth = .2;%.16
fG = -3;
plotOps = {'wing2fot_histogram';'wing2fot_histogram_smoothed' % 1 - 2
    'escape';'spaghettiPlot_zeroFly';'departureScatter' % 3 - 5
    'departureHistogram_zeroFly';'departureHistogram_zeroStim'}; % 6 - 7
graphName = plotOps{plotChoice};
dirOpsA = {'Jump','Departure','Trajectory','Meatball'};
dirOpsB = {'zeroStimAzi','zeroFlyAzi','3D','elevation'};
dirOpsC = {'Angle','Dist','Binomial'};
dirOpsD = {'RhoVar','RhoOne'};
segmentOps = {'FullJumpDur','Fot','MidLegExt','FullLegExt','FullDur','MidLoom','JumpCap','LateLoom','MvntCap'};
if strcmp(graphName,'spaghettiPlot_zeroFly')
    graphName = cat(2,graphName,'_',segmentOps{segmentChoice},'_overlay');
elseif strcmp(graphName,'escape')
    graphName = cat(2,graphName,dirOpsA{dirChoice(1)},'_',dirOpsB{dirChoice(2)},...
        '_',dirOpsC{dirChoice(3)},dirOpsD{dirChoice(4)});
end
optionsPath = 'Z:\Data_pez3000_analyzed\WRW_graphing_variables\graphOptions.mat';
%%%%%
makeGraphOptionsStruct([],graphName)
[plotTable,excelTable] = addPlotData(graphTable);
%%%%%
if numel(get(0,'children')) > 10
    close all
end

graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
sheetName2plot = graphOptions.sheetName;
graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
groupIDlist = graphOptions.groupIDlist;
plotIDlist = graphOptions.plotIDlist;
xCt = numel(plotIDlist);
grpID = groupIDlist{groupChoice};
grpBool = strcmp(grpID,plotTable.groupID);
if max(strcmp(excelTable.Properties.VariableNames,'repeat'))
    grpBool = grpBool | plotTable.repeat;
end
grpTable = plotTable(grpBool,:);
[~,grpOrder] = sort(grpTable.order);
grpTable = grpTable(grpOrder,:);
grpIDlist = grpTable.Properties.RowNames;
dataRef = find(ismember((1:xCt),grpTable.order));
subset2plot = graphOptions.subset2plot;
xCt = numel(dataRef);
% return

if xCt > 33
    xCt = 33;
    disp('overcount!')
end

%group colors
colorOps = {[55 126 184];[228 26 28];[77 175 74];[152 78 163];[255 127 0]};
colorOpsA = cellfun(@(x) (x/255).^(0.7),colorOps,'uniformoutput',false);
colorOps = {[55 126 184];[228 26 28];[77 175 74];[152 78 163];[255 127 0]};
colorOpsB = cellfun(@(x) (x/255).^(1.2),colorOps,'uniformoutput',false);
fontSa = 18; fontSb = 12; fontSc = 12;
fontSa = fontSa+fG; fontSb = fontSb+fG; fontSc = fontSc+fG;
colorType = 1;
if colorType == 1
    fontC = [0 0 0];
    backC = [1 1 1];
    lineCa = [0 0 0];
    lineCb = [0.5 0.5 0.5]+0.2;
    colorOps = cat(1,colorOpsB,colorOpsA);
else
    fontC = [1 1 1];
    backC = [0 0 0];
    lineCa = [0 0 0];
    lineCb = [0 0 0]+0.5;
    colorOps = cat(1,colorOpsA,colorOpsB);
end

figure
figPos = [1988 43 1778 1030];
set(gcf,'position',figPos)
hax = zeros(1,xCt);
hax(1) = axes;
haxPos = get(hax(1),'position');
growShiftAx = [0.1,-0.1,0,0.02];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
if ~isempty(strfind(graphName,'wing2fot_histogram'))
    graphName = [graphName '_gaussFitCt-' num2str(gaussFitCt) '_barWidth-' regexprep(num2str(barwidth),'\.','pt')];
    if xCt > 15
        rows = 4;
        axCt = xCt;
    elseif xCt > 6
        rows = 3;
        axCt = 15;
    else
        rows = 2;
        axCt = 6;
    end
    
    if axCt < 4, axCt = 4; end
    axWH = [haxPos(3)/ceil(axCt/rows) haxPos(4)/rows];
    startPosY = repmat(fliplr(axWH(2)*(0:rows-1)),ceil(axCt/rows),1);
    startPosX = repmat(axWH(1)*(0:ceil(axCt/rows)-1),1,rows);
    options = statset('MaxIter',500);
    xmin = 0.0;      xmax = 7;
    intervals = (xmin:barwidth:xmax);
    fitX = linspace(xmin,xmax,1000);
%     xtickval = unique(round(logspace(log10(1),log10(300),10)/5)*5);
    xtickval = unique(round(logspace(log10(1),log10(300),10)/10)*10)
    xtickval = [1.25 2.5 5 xtickval(2:end)];
    xtickpos = log(xtickval);
    
%     xmin = 0;      xmax = 300;
%     intervals = (xmin:barwidth:xmax);
%     fitX = linspace(xmin,xmax,1000);
%     xtickval = round(xmin:50:xmax);
%     xtickpos = xtickval;
    ytickval = [0 1]; ytickpos = [0 1];
    xAxisName = 'ms';
    grpTable.wing2fot_hist = cell(size(grpTable,1),1);
elseif ~isempty(strfind(graphName,'spaghettiPlot')) || ~isempty(strfind(graphName,'escape'))
    if xCt > 15
        rows = 4;
    elseif xCt > 8
        rows = 3;
    else
        rows = 2;
    end
    axWH = [haxPos(3)/ceil(xCt/rows) haxPos(4)/rows];
    startPosY = repmat(fliplr(axWH(2)*(0:rows-1)),ceil(xCt/rows),1);
    startPosX = repmat(axWH(1)*(0:ceil(xCt/rows)-1),1,rows);
    xAxisName = 'mm';
    if ~isempty(strfind(graphName,'escape'))
        if ~isempty(strfind(graphName,'Trajectory'))
            %         xtickpos = (-3:1.5:3);
            xtickpos = linspace(-1.5,1.5,3);
        elseif ~isempty(strfind(graphName,'Departure'))
            %         xtickpos = (-0.5:0.25:0.5);
            xtickpos = linspace(-1,1,5);
        else
            xtickpos = linspace(-1,1,3);
        end
    else
        xtickpos = (-0.5:0.5:0.5);
    end
%     xtickpos = xtickpos/2;
    xtickval = xtickpos;
    if ~isempty(strfind(graphName,'elevation'))
        ytickpos = xtickpos(xtickpos >= 0);
    else
        ytickpos = xtickpos;
    end
    ytickval = ytickpos;
elseif ~isempty(strfind(graphName,'departureScatter'))
%     xtickpos = (-180:60:180);
    xtickpos = (0:30:180);
    xtickval = xtickpos;
    ytickpos = linspace(-360,180,7);
    ytickval = ytickpos+360;
    axWH = [haxPos(3) haxPos(4)];
    startPosY = axWH(2)*0.1;
    startPosX = axWH(1)*0.1;
    if ~isempty(strfind(graphName,'zeroFly'))
    elseif ~isempty(strfind(graphName,'zeroStim'))
    end
    xAxisName = 'Stimulus Azimuth (degrees)';
end
%%%% %Title and other labels
if strcmp(groupIDlist{groupChoice},sheetName2plot)
    pdfName = [sheetName2plot '_' graphName '_' graphOptions.subset2plot];
else
    pdfName = [sheetName2plot '_' graphName '_' groupIDlist{groupChoice} '_' graphOptions.subset2plot];
end
htit = text(xtickpos(1),ytickpos(end)*1.3,pdfName,...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',round(fontSa*1.2));

if ~isempty(strfind(graphName,'departureScatter'))
    singleAx = true;
    axWH = [haxPos(3) haxPos(4)];
    startPosY = axWH(2)*0.1;
    startPosX = axWH(1)*0.1;
end
thetaCell = cell(xCt,1);
alphaCell = cell(xCt,1);
h4leg = zeros(xCt,1);

for iterE = 1:xCt
%     if ~dataRef(iterE)
%         if iterE == 1
%             delete(hax(1))
%         end
%         continue
%     end
    ndxName = grpIDlist{iterE};
    data = grpTable.returnData(ndxName)';
    try
        data = (cat(1,data{:}))';
    catch
        data = (cat(2,data{:}))';
    end
    if iterE > 1
        if ~singleAx
            hax(iterE) = axes;
        else
            hax(iterE) = hax(iterE-1);
            hold on
        end
    end
    if hax(iterE) == 0
        continue
    end
    set(hax(iterE),'nextplot','add')
    if isempty(data)
%         if singleAx
            delete(hax(iterE))
%         end
        continue
    end
    if ~isempty(strfind(graphName,'wing2fot_histogram'))
        data(data > xmax) = [];
        ctMax = Inf;
        if numel(data) > ctMax
            randVec = randperm(ctMax);
            data = data(randVec(1:ctMax));
        end
        [f,x] = hist(data,intervals);
        if ~isempty(strfind(graphName,'smoothed'))
            f = f/range(f);
            f = round(smooth(f,3)*100);
            fitdata = cell(numel(f),1);
            for iterF = 1:numel(x)
                fitdata{iterF} = repmat(x(iterF),f(iterF),1);
            end
            fitdata = cat(1,fitdata{:});
        else
            fitdata = data(:);
        end
        fitdata(isnan(fitdata)) = [];
        fitdata(fitdata < 1) = [];
        jumpPct = grpTable.jumpCount(ndxName)/grpTable.plotCount(ndxName);
        f = f/range(f);
%         f = f*(jumpPct^0.15);
        bar(x,f,'stack','barwidth',1,'facecolor',colorOps{1}.^0.7);
        orig_state = warning;
        warning('off','all')
        
        obj = cell(numel(gaussFitCt),1);
        aic = zeros(numel(gaussFitCt),1);
        for iters = 1:numel(gaussFitCt)
            conv = false;
            iterations = gaussFitCt(iters)+1;
            while ~conv && iterations > 1
                iterations = iterations - 1;     int_cond = [];
                for iterI = 1:(iterations - 1);
                    int_cond = [int_cond;ones(round(length(fitdata)/iterations),1)*iterI];
                end
                if isempty(iterI), iterI = 0; end
                int_cond = [int_cond;ones((length(fitdata) - length(int_cond)),1)*(iterI+1)];
                try
                    %                 obj = gmdistribution.fit(fitdata,iterations,'Options',options,'Start',int_cond);
                    obj{iters} = fitgmdist(fitdata,iterations,'Options',options,'Start',int_cond);
                    conv = obj{iters}.Converged;
                    aic(iters) = obj{iters}.BIC;
%                     obj.BIC
%                     obj.NegativeLogLikelihood
                catch ME
%                                     getReport(ME)
                end
            end
        end
        warning(orig_state)
        obj = obj{aic == min(aic)};
        if ~isempty(obj) && obj.Converged
            fitY = pdf(obj,fitX');
            fitY = fitY/range(fitY);
            plot(fitX,fitY,'color',colorOps{2}.^1.3,'LineWidth',3);
        end
        plot(log([41/6 41/6]),[0 1],'color','k','linewidth',3,'linestyle','--');
        ycolor = 'none';
        grpTable.wing2fot_hist(ndxName) = {[x(:),f(:)]};
    elseif ~isempty(strfind(graphName,'spaghettiPlot'))
        axis equal
        plot([xtickpos(1) xtickpos(end)],zeros(1,2),'color',lineCb)
        plot(repmat(median(xtickpos),1,2),[ytickpos(1) ytickpos(end)],'color',lineCb)
        comboPlot = 3;
        brkPts = find(isnan(data(1,:)));
        brkCt = numel(brkPts);
        if comboPlot == 1
            if brkCt > 100, brkPt = brkPts(100); else brkPt = brkPts(end); end
            plot(data(1,1:brkPt),data(2,1:brkPt),'color',colorOps{3},'LineWidth',0.5);
        elseif comboPlot == 2
            if brkCt > 100, brkCt = 100; end
            for iterBrk = 1:brkCt
                if iterBrk == 1
                    data2plot = data(:,1:brkPts(iterBrk));
                else
                    data2plot = data(:,brkPts(iterBrk-1)+1:brkPts(iterBrk));
                end
                plot(data2plot(1,:),data2plot(2,:),'LineWidth',0.5)%,'color',cA{iterC});
            end
        else
            if brkCt > 100, brkCt = 100; end
            data2plot = [];
            for iterBrk = 1:brkCt
                if grpTable.jumpTest{ndxName}(iterBrk) == 1
                    continue
                end
                if iterBrk == 1
                    data2plot = cat(2,data2plot,data(:,1:brkPts(iterBrk)));
                else
                    data2plot = cat(2,data2plot,data(:,brkPts(iterBrk-1)+1:brkPts(iterBrk)));
                end
            end
            if ~isempty(data2plot)
                hplot = plot(data2plot(1,:),data2plot(2,:),'LineWidth',1,...
                    'linestyle','-','color',colorOps{3});
                hplot.Color(4) = 0.5;
            end
            data2plot = [];
            for iterBrk = 1:brkCt
                if grpTable.jumpTest{ndxName}(iterBrk) == 0
                    continue
                end
                if iterBrk == 1
                    data2plot = cat(2,data2plot,data(:,1:brkPts(iterBrk)));
                else
                    data2plot = cat(2,data2plot,data(:,brkPts(iterBrk-1)+1:brkPts(iterBrk)));
                end
            end
            if ~isempty(data2plot)
                hplot = plot(data2plot(1,:),data2plot(2,:),'LineWidth',1,...
                    'linestyle','-','color',colorOps{4});
                hplot.Color(4) = 0.5;
            end
        end
        [stimx,stimy] = pol2cart(grpTable.azimuth(ndxName)*(pi/180),max(ytickpos)/2);
        if ~strcmp(ndxName,'noStim') && ~strcmp(grpTable.stimType{ndxName},'Photoactivation')
            stimX = [0 stimx];
            stimY = [0 stimy];
        else
            stimX = [0 0];
            stimY = [NaN NaN];
        end
        plot(stimX,stimY,'color',colorOps{1},'LineWidth',2);
        ycolor = fontC;
    elseif ~isempty(strfind(graphName,'escape'))
        axis equal
        plot([xtickpos(1) xtickpos(end)],zeros(1,2),'color',lineCb,'linewidth',2)
        plot(repmat(median(xtickpos),1,2),[ytickpos(1) ytickpos(end)],'color',lineCb,'linewidth',2)
        [midX,midY] = pol2cart(pi/4,ytickpos(end));
        if ~isempty(strfind(graphName,'elevation'))
            midY = cat(1,midY,0,midY);
            midX = cat(1,-midX,0,midX);
        else
            midY = cat(1,midY,-midY,NaN,-midY,midY);
            midX = cat(1,-midX,midX,NaN,-midX,midX);
        end
        plot(midX,midY,'color',lineCb,'linewidth',2)
        plot(-midX,midY,'color',lineCb,'linewidth',2)
        if ~isempty(strfind(graphName,'elevation'))
            circT = linspace(0,pi,500);
        else
            circT = linspace(-pi,pi,500);
        end
        [circX,circY] = pol2cart([circT circT],[zeros(size(circT))+1 zeros(size(circT))+0.5]);
        plot(circX,circY,'color',lineCb,'LineWidth',2,'linestyle','-');
        if graphOptions.sliceOption > 1
            sliceOption = graphOptions.sliceOption;
            aziGrpCt = sliceOption;% 5 or 7 is recommended ... 13 works too
            if aziGrpCt/2 == round(aziGrpCt/2)
                aziBreaks = linspace(0,180,aziGrpCt*2+1);
                aziNames = aziBreaks(2:2:aziGrpCt*2+1);
                aziBreaks = aziBreaks(1:2:aziGrpCt*2+1);
            else
                aziBreaks = linspace(0,180,aziGrpCt*2-1);
                aziNames = aziBreaks(1:2:aziGrpCt*2-1);
                aziBreaks = [aziBreaks(1) aziBreaks(1,2:2:aziGrpCt*2-1) aziBreaks(end)];
            end
            brkRef = find(aziNames == grpTable.azimuth(ndxName));
            stimT = aziBreaks(brkRef:brkRef+1)*(pi/180);
            stimT = linspace(stimT(1),stimT(2),100);
            [stimX,stimY] = pol2cart(stimT,zeros(size(stimT))+max(ytickpos)*0.75);
            plot(stimX,stimY,'color',[0 0 0]+0.5,'linewidth',5);
        end
        
        dirData = data(1,:);
        eleData = -data(3,:);
        if strcmp(grpTable.stimType{ndxName},'Photoactivation')
%             dirData = abs(dirData);
        end
%         rhoVec = data(2,:);
%         if ~isempty(strfind(graphName,'RhoOne'))
%             dirData(rhoVec < rhoThresh) = [];
%             grpTable.jumpTest{ndxName}(rhoVec < rhoThresh) = [];
%             rhoVec(rhoVec < rhoThresh) = [];
%             rhoVec = rhoVec*0+1;
%         end
        
        rhoVec = data(4,:);
        dirData(rhoVec < rhoThresh) = [];
        eleData(rhoVec < rhoThresh) = [];
        grpTable.jumpTest{ndxName}(rhoVec < rhoThresh) = [];
        rhoVec(rhoVec < rhoThresh) = [];
        rhoVec = rhoVec*0+1;
        [dataxTot,datayTot,datazTot] = sph2cart(deg2rad(dirData),deg2rad(-eleData),rhoVec);
        
%         [dataxTot,datayTot] = pol2cart(dirData*(pi/180),rhoVec);
        for iterJ = 1:2
            if iterJ == 1
                datax = dataxTot(grpTable.jumpTest{ndxName});
                if ~isempty(strfind(graphName,'elevation'))
                    datay = datazTot(grpTable.jumpTest{ndxName});
                else
                    datay = datayTot(grpTable.jumpTest{ndxName});
                end
                colorRef = 3;
            else
                if strcmp(subset2plot,'jumping')
                    continue
                end
                datax = dataxTot(~grpTable.jumpTest{ndxName});
                if ~isempty(strfind(graphName,'elevation'))
                    datay = datazTot(~grpTable.jumpTest{ndxName});
                else
                    datay = datayTot(~grpTable.jumpTest{ndxName});
                end
                colorRef = 4;
            end
            datax(isnan(datax)) = [];
            datay(isnan(datay)) = [];
            if numel(datax > 5)
                if plotCircArrow == 1
                    avgx = median(datax(:));
                    avgy = median(datay(:));
                    [avgT,avgR] = cart2pol(avgx,avgy);
                    if avgR >= 0.3
                        avgT = cat(2,avgT,avgT+pi/45*[1 -1],avgT);
                        avgR = max(ytickpos);
                        avgR = cat(2,avgR,[avgR avgR]*(0.7),avgR);
                        [avgx,avgy] = pol2cart(avgT,avgR);
                        plot(avgx,avgy,'color',colorOps{colorRef},'linewidth',2);
                    end
                end
                if showDotsOrLines == 1
                    plotResample = randperm(numel(datay));
                    ctMax = 150;
                    if numel(plotResample) < ctMax
                        ctMax = numel(plotResample);
                    end
                    for iterP = 1:ctMax
                        if ~isempty(strfind(graphName,'elevation'))
                            plot([0 datax(plotResample(iterP))],[0 datay(plotResample(iterP))],...
                                'color',colorOps{colorRef},'linewidth',1)
                        else
                            plot(datax(plotResample(iterP)),datay(plotResample(iterP)),...
                                '.','color',colorOps{colorRef},'markersize',6)
                        end
                    end
                end
                if plotCircHist == 1
                    
                    if iterJ == 1
                        histData = dirData(grpTable.jumpTest{ndxName});
                    else
                        histData = dirData(~grpTable.jumpTest{ndxName});
                    end
                    if ~isempty(strfind(graphName,'elevation'))
                        aziGrpCt = 60;% 9 or 13 or 25 or 31
                    else
%                         aziGrpCt = 25;% 9 or 13 or 25 or 31
                        aziGrpCt = 13;% 9 or 13 or 18 or 25 or 31
                    end
                    aziBreaks = linspace(-180,180,aziGrpCt*2-1);
                    aziNames = aziBreaks(1:2:aziGrpCt*2-1);
                    aziBreaks = [aziBreaks(1) aziBreaks(1,2:2:aziGrpCt*2-1) aziBreaks(end)];
                    dataCt = numel(histData);
                    histVec = zeros(numel(aziBreaks)-1,1);
                    for iterP = 1:numel(aziBreaks)-1
                        histBool = histData > aziBreaks(iterP) & histData < aziBreaks(iterP+1);
                        if iterP == 1
                            histBoolB = histData > aziBreaks(end-1) & histData < aziBreaks(end);
                            histBool = histBool | histBoolB;
                        end
                        histVec(iterP) = sum(histBool)/dataCt;
                    end
                    histVec = histVec/max(histVec);
                    histVec = cat(2,histVec(:),histVec(:))';
                    aziBreaks = cat(2,aziBreaks(:),aziBreaks(:))';
                    aziBreaks = circshift(aziBreaks(:),[1 0]);
                    if iterJ == 1
                        aziBreaks = aziBreaks+2;
                        histVec = histVec*0.92;
                    else
                        aziBreaks = aziBreaks-2;
                        histVec = histVec*0.88;
                    end
                    [plotX,plotY] = pol2cart(deg2rad(aziBreaks(3:end)),histVec(:));
                    plotX = cat(1,plotX(:),plotX(1));
                    plotY = cat(1,plotY(:),plotY(1));
                    plot(plotX,plotY,'color',colorOps{colorRef},'linewidth',2)
                end
                if ~isempty(strfind(graphName,'elevation'))
                    medData = cart2pol(datax,datay);
                    q1 = prctile(medData,25); q2 = prctile(medData,50);
                    q3 = prctile(medData,75); n = numel(medData);
                    err = 1.57*(q3-q1)/sqrt(n);
                    errR = max(ytickpos);
                    confT = linspace(q2-err,q2+err,100);
                    [errX,errY] = pol2cart([confT NaN q2 q2],[zeros(size(confT))+errR*1.025 NaN errR*1.0 errR*1.05]);
                    plot(errX,errY,'k','linewidth',2)
                else
                    medVals = zeros(2,2);
                    for iterErr = 1:2
                        if iterErr == 1
                            medData = datax;
                        else
                            medData = datay;
                        end
                        q1 = prctile(medData,25); q2 = prctile(medData,50);
                        q3 = prctile(medData,75); n = numel(medData);
                        err = 1.57*(q3-q1)/sqrt(n);
                        medVals(iterErr,1) = q2;
                        medVals(iterErr,2) = err;
                    end
                    if circConfStyle == 1
                        plot([medVals(1,1)-medVals(1,2) medVals(1,1)+medVals(1,2)],...
                            [medVals(2,1) medVals(2,1)],'linewidth',1,'color','k')
                        plot([medVals(1,1) medVals(1,1)],...
                            [medVals(2,1)-medVals(2,2) medVals(2,1)+medVals(2,2)],...
                            'linewidth',1,'color','k')
                    elseif circConfStyle == 2
                        xCenter = medVals(1,1);
                        yCenter = medVals(2,1);
                        xRadius = medVals(1,2);
                        yRadius = medVals(2,2);
                        theta = 0 : 0.01 : 2*pi;
                        x = xRadius * cos(theta) + xCenter;
                        y = yRadius * sin(theta) + yCenter;
                        fill(x,y,'k','LineWidth',2);
                    end
                end
%                 if isempty(strfind(graphName,'Meatball'))
%                     plot([0 avgx],[0 avgy],'color','k','LineWidth',5);
%                     plot([0 avgx],[0 avgy],'color',colorOps{2},'LineWidth',3);
%                 end
            end
        end
        
        ycolor = fontC;
        if graphOptions.sliceOption == 1
            [stimx,stimy] = pol2cart(grpTable.azimuth(ndxName)*(pi/180),max(ytickpos));
            if ~strcmp(ndxName,'noStim') && ~strcmp(grpTable.stimType{ndxName},'Photoactivation')
                stimX = [0 stimx];
                stimY = [0 stimy];
            else
                stimX = [0 0];
                stimY = [NaN NaN];
            end
            plot(stimX,stimY,'color',colorOps{1},'LineWidth',3);
        end
        
    elseif ~isempty(strfind(graphName,'departureScatter'))
%         xbaseX = [xtickpos;xtickpos;xtickpos];
%         xbaseY = [zeros(size(xtickpos))+ytickpos(1);zeros(size(xtickpos))+ytickpos(end);zeros(size(xtickpos))+ytickpos(1)];
%         plot(xbaseX(:),xbaseY(:),'linewidth',1,'color',lineCb)
%         ybaseY = [ytickpos;ytickpos;ytickpos];
%         ybaseX = [zeros(size(ytickpos))+xtickpos(1);zeros(size(ytickpos))+xtickpos(end);zeros(size(ytickpos))+xtickpos(1)];
%         plot(ybaseX(:),ybaseY(:),'linewidth',1,'color',lineCb)
        
        theta = data(2,:);
        alpha = data(1,:);
        plot(theta,alpha,'.','color',colorOps{3},'markersize',12)
        ycolor = fontC;
        thetaCell{iterE} = theta;
        alphaCell{iterE} = alpha;
%         stimx = max(ytickpos)/2;
%         stimy = 0;
%         if ~isempty(strfind(graphName,'zeroFly'))
%             [stimx,stimy] = pol2cart(str2double(labelSpl{2})*(pi/180),stimx);
%         end
%         plot([0 stimx],[0 stimy],'color',cA{2},'LineWidth',2);
    end
    
    if iterE == 1 || ~singleAx
        haxPos = get(hax(iterE),'position');
        haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
            haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
        haxPos = [haxPos(1:2)+[startPosX(iterE) startPosY(iterE)] axWH.*[0.85 0.75]];
        
        set(hax(iterE),'position',haxPos,'color','none','box','off')
        %     continue
        set(hax(iterE),'xticklabel',[],'position',haxPos,'color','none','box','off',...
            'Xlim',[xtickpos(1)-range(xtickpos)*0.1 xtickpos(end)+range(xtickpos)*0.1],...
            'YLim',[ytickpos(1)-range(ytickpos)*0.1 ytickpos(end)+range(ytickpos)*0.1],...
            'tickdir','out','xtick',[],'yticklabel',[],'ytick',ytickpos,'ycolor','none','xcolor','none')
        if ~isempty(strfind(graphName,'escape'))
            aziVals = linspace(0,180,5);
            rhoVals = zeros(size(aziVals))+max(xtickpos);
            [xpos,ypos] = pol2cart(deg2rad(aziVals),rhoVals*1.1);
            for iterL = 1:numel(ypos)
                text(xpos(iterL),ypos(iterL),num2str(aziVals(iterL)),'rotation',aziVals(iterL)-90,...
                    'color',ycolor,'horizontalalignment','center','fontsize',fontSb);
            end
        else
            for iterL = 1:numel(ytickval)
                text(xtickpos(1)-range(xtickpos)*0.05,ytickpos(iterL),num2str(ytickval(iterL)),'rotation',0,...
                    'color',ycolor,'horizontalalignment','right','fontsize',fontSb);
            end
        end
        xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.2];
        if isempty(strfind(graphName,'departureScatter'))
            htxtA = text(xtickpos(1),xlabelPos(2),strtrim(grpTable.plotID{ndxName}),...
                'horizontalalignment','left','rotation',0,...
                'color',fontC,'interpreter','none','fontsize',fontSc);
            jumpCt = sum(grpTable.jumpTest{ndxName});
            plotCt = grpTable.plotCount(ndxName);
            if strcmp(subset2plot,'jumping')
                dataCountTotal = ['n = ' num2str(jumpCt)];
            elseif strcmp(subset2plot,'nojump')
                dataCountTotal = ['n = ' num2str(plotCt)];
            else
                dataCountTotal = ['n=' num2str(plotCt),'   ' num2str(jumpCt),...
                    ' jump   ' num2str(plotCt-jumpCt) ' nojump'];
            end
        else
            dataCountTotal = ['n = ' num2str(sum(grpTable.jumpCount(:)))];
            ylabelPos = [xtickpos(1)-range(xtickpos)*0.2 mean(ytickpos) 1];
            htxtA = text(ylabelPos(1),ylabelPos(2),'Departure Azimuth (degrees)',...
                'horizontalalignment','center','rotation',90,...
                'color',fontC,'interpreter','none','fontsize',fontSa);
        end
        xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.25];
        htxtB = text(xtickpos(1),xlabelPos(2),dataCountTotal,...
            'horizontalalignment','left','rotation',0,'VerticalAlignment','top',...
            'color',fontC,'interpreter','none','fontsize',fontSb,'parent',hax(iterE));
    end
    if iterE == xCt && isempty(strfind(graphName,'escape'))
        htxtA.Position(2) = htxtA.Position(2)-range(ytickpos)*0.25;
        htxtB.Position(2) = htxtB.Position(2)-range(ytickpos)*0.25;
        set(get(hax(iterE),'xlabel'),'units','normalized')
        xlabelPos = get(get(hax(iterE),'xlabel'),'position');
        xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.2];
        for iterL = 1:numel(xtickval)
            if isnan(xtickval(iterL)), continue, end
            text(xtickpos(iterL),xlabelPos(2),num2str(xtickval(iterL)),'rotation',0,...
                'color',fontC,'horizontalalignment','center','fontsize',fontSb);
        end
        xlabelPos = get(get(hax(iterE),'xlabel'),'position');
        xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.3];
        text(xlabelPos(1),xlabelPos(2),xAxisName,...
            'horizontalalignment','center','rotation',0,...
            'color',fontC,'interpreter','none','fontsize',fontSb);
    end
end
scatterFit = 3;
if ~isempty(strfind(graphName,'departureScatter'))
    theta = cat(2,thetaCell{:});
    alpha = cat(2,alphaCell{:});
    if scatterFit == 1
        fitLeg = 3;
    else
        fitLeg = 3;
    end
    fitPts = zeros((max(xtickpos)-fitLeg*2),3);
    fitX = zeros((max(xtickpos)-fitLeg*2),1);
    fitZ = zeros((max(xtickpos)-fitLeg*2),1);
    for iterF = 1:(max(xtickpos)-fitLeg*2)
        binStop = iterF+fitLeg*2+1;
        thetaNdcs = theta > iterF & theta <= binStop;
        q2 = median(alpha(thetaNdcs));
        q1 = prctile(alpha(thetaNdcs),25);
        q3 = prctile(alpha(thetaNdcs),75);
        fitPts(iterF,:) = [q1 q2 q3];
        fitX(iterF) = iterF+fitLeg+1;
        fitZ(iterF) = median(theta(thetaNdcs));
    end
    if scatterFit == 2 || scatterFit == 3
        for iterF = 2%1:3
            plot(fitX,fitPts(:,iterF),'color',[0 0 0.5],'linewidth',3)
        end
    end
    if scatterFit == 1 || scatterFit == 3
        for iterF = 2%1:3
            fitY = fitPts(:,iterF);
            cfun = fit(fitX,fitY,'power2');
            plotX = (min(xtickpos):max(xtickpos));
            plot(plotX,feval(cfun,plotX),'color',[0.5 0 0],'linewidth',3)
        end
    end
    plot([0 180],[-180 0],'color',fontC,'linewidth',1,'linestyle','-')
    plot([0 180],[-360 -180],'color',fontC,'linewidth',1,'linestyle','-')
    plot([0 180],[0 180],'color',fontC,'linewidth',1,'linestyle','-')
    plot([0 180],[-180 -180],'color',fontC,'linewidth',2,'linestyle','--')
    plot([0 180],[0 0],'color',fontC,'linewidth',2,'linestyle','--')
elseif singleAx == 1
    hLeg = legend(h4leg,grpIDlist,'interpreter','none');
    legPos = get(hLeg,'pos');
    legPos(1) = legPos(1)*1.3;
    set(hLeg,'pos',legPos)
end
exportVisualizationFigure(pdfName,saveFig)
