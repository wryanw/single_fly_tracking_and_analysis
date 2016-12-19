function [pdfName,plotTable] = graphing_pez3000_vectorPlots(graphTable)
%%
if ~exist('graphTable','var') || isempty(mfilename)
    graphTable = makeGraphingTable;
end
%%

clearvars -except graphTable
rng('default');
rng(19);

% use with Motion
segmentChoice = 4; % 1-full 10 ms...2-fot...3-mid leg extend...4-full leg
    % extend...5-entire vector...6-half loom...7-jumpCap...8-last half
    % loom...9-distCap...10-fullStim...11-maxVelWin
unitChoice = 4; % 1 - velocity ; 2 - acceleration ; 3 - distance ; 4 - force
% use with relMotion
motionChoice = 2; % 1-movement 2-walking 3-turning 4-sidepass

plotChoice = 1;
grpChoice = 1; % if greater than '0', grp to use
smoothData = 0;
curveFit = [0 4];
saveFig = 2; % 1-pdf ; 2-eps ; 3-png
singleAx = 0;
medianAndBounds = 1;% 1 - show median trace ; 2 - show iqr traces
leaveSpaces4omits = 0;
showCountTrace = 1;
xScale = 0; % 0-normal ; 1-log
showMedian = 1; % for FOT plots
barWidth_ms = 27; % 27 was original
errorMode = [2 1];
errorBinSize = 100;
explicitX = [1 0 4 5];
% explicitX = [0 -500 100 7];
% explicitX = [1 0 500 3];
padDir = -1; % 1 - pad back end ; -1 - pad front end
% explicitY = [1 -60 60 3];
explicitY = [1 0 250 6];
% explicitY = [1 0 1.5 4];
% explicitY = [1 0 15 3];
fG = -3; % font growth
explicitTimeBounds = [0 1 500]; % 1 to use, start time, stop time
plotOps = {'absMotion';'relMotion';'binaryMovement' % 1 - 3
    'frameOfX_time_TOsum';'frameOfX_time_TOhist';'frameOfX_time_MVNTsum' % 4 - 6
    'frameOfX_time_TOraster';'jaabaReaching'}; % 7 - 8
graphName = plotOps{plotChoice};
motionOps = {'movement','walking','turning','sidepass'};
segmentOps = {'FullJumpDur','Fot','MidLegExt','FullLegExt','FullDur','MidLoom',...
    'JumpCap','LateLoomMvntCap','MvntCap','FullStim','MaxVelWin'};
unitOps = {'Vel','Accel','Dist','Force'};
if strcmp(graphName,'absMotion') && motionChoice ~= 3
    motionChoice = 1;
end
if strcmp(graphName,'jaabaReaching')
    showCountTrace = 0;
    smoothData = 1;
    explicitY = [1 0 1 3];
end
binaryMovement = 0;
if explicitTimeBounds(1) == 1
    segmentOps{segmentChoice} = ['timeBounds' num2str(explicitTimeBounds(2)) 'to' num2str(explicitTimeBounds(3))];
end
if strcmp(graphName,'binaryMovement')
    graphName = cat(2,'absMotion_',motionOps{1},'_',unitOps{1},'_',segmentOps{segmentChoice},'_vector');
    binaryMovement = 1;
    smoothData = 1;
    medianAndBounds = 0;
    showCountTrace = 1;
    explicitY(1) = 0;
elseif ~isempty(strfind(graphName,'frameOfX'))
    smoothData = 0;
    medianAndBounds = 0;
    showCountTrace = 0;
    explicitY(1) = 0;
else
    graphName = cat(2,graphName,'_',motionOps{motionChoice},'_',unitOps{unitChoice},'_',segmentOps{segmentChoice},'_vector');
end

optionsPath = 'Z:\Data_pez3000_analyzed\WRW_graphing_variables\graphOptions.mat';
%%%%%
makeGraphOptionsStruct([],graphName)
graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
graphOptions.padOption = padDir;
save(optionsPath,'graphOptions')
[plotTable,excelTable] = addPlotData(graphTable);
%%%%%
graphOptionsLoading = load(optionsPath);
graphOptions = graphOptionsLoading.graphOptions;
sheetName2plot = graphOptions.sheetName;
groupIDlist = graphOptions.groupIDlist;
grpCt = numel(groupIDlist);
if grpChoice > 0
    grpOps = grpChoice;
else
    close all
    grpOps = (1:grpCt);
end
plotIDlist = graphOptions.plotIDlist;

if medianAndBounds == 1
    graphName = cat(2,graphName,'_median');
elseif medianAndBounds == 2
    graphName = cat(2,graphName,'_iqrOnly');
end
if numel(get(0,'children')) > 10
    close all
end
if binaryMovement == 1
    graphName = regexprep(graphName,'absMotion','binaryMovement');
end
for iterG = grpOps
    grpID = groupIDlist{iterG};
    grpBool = strcmp(grpID,plotTable.groupID);
    if max(strcmp(excelTable.Properties.VariableNames,'repeat'))
        grpBool = grpBool | plotTable.repeat;
    end
    grpTable = plotTable(grpBool,:);
    [~,grpOrder] = sort(grpTable.order);
    grpTable = grpTable(grpOrder,:);
    grpIDlist = grpTable.Properties.RowNames;
    if leaveSpaces4omits == 1
        xCt = numel(plotIDlist);
        dataRef = ismember((1:xCt),grpTable.order);
    else
        xCt = numel(grpIDlist);
        dataRef = grpTable.order;
    end
    subset2plot = graphOptions.subset2plot;
    
    % return
    
    if xCt > 40
        xCt = 40;
        disp('overcount!')
    end
    
    %group colors
    colorOps = {[55 126 184];[228 26 28];[77 175 74];[152 78 163];[255 127 0]};
    colorOpsA = cellfun(@(x) (x/255).^(0.7),colorOps,'uniformoutput',false);
    colorOps = {[55 126 184];[228 26 28];[77 175 74];[152 78 163];[255 127 0]};
    colorOpsB = cellfun(@(x) (x/255).^(1.2),colorOps,'uniformoutput',false);
    fontSa = 16; fontSb = 12; fontSc = 12;
    fontSa = fontSa+fG; fontSb = fontSb+fG; fontSc = fontSc+fG;
    colorType = 1;
    if colorType == 1
        fontC = [0 0 0];
        backC = [1 1 1];
        lineCb = [0.5 0.5 0.5]+0.2;
        colorOps = cat(1,colorOpsB,colorOpsA);
    else
        fontC = [1 1 1];
        backC = [0 0 0];
        lineCb = [0 0 0]+0.5;
        colorOps = cat(1,colorOpsA,colorOpsB);
    end
    
    figure
    if size(get(0,'monitorpositions'),1) == 2
        figPos = [1988 43 1778 1030];
    else
        figPos = [88 43 1778 1030];
    end
    set(gcf,'color',backC,'position',figPos)
    hax = zeros(1,xCt);
    hax(1) = axes;
    haxPos = get(hax(1),'position');
    growShiftAx = [0.1,-0.1,0,0.02];%grow width, grow height, shift X, shift Y
    haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
        haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
    
    %%%% %Title and other labels
    pdfName = [sheetName2plot '_' graphName '_' groupIDlist{iterG} '_' graphOptions.subset2plot];
    if singleAx == 1
        cols = 2;
        rows = 2;
    else
        cols = ceil(xCt/6);
        rows = ceil(xCt/cols);
    end
    if cols == 1
        figPos(3) = figPos(3)*0.4;
        figPos(1) = figPos(1)+200;
        set(gcf,'pos',figPos)
    elseif xCt <= 12
        figPos(3) = figPos(3)*0.7;
        figPos(4) = figPos(4)*0.8;
        figPos(1) = figPos(1)+300;
        set(gcf,'pos',figPos)
    end
    axWH = [haxPos(3)/cols haxPos(4)/rows];
    startPosY = repmat(fliplr(axWH(2)*(0:rows-1)),cols,1)';
    startPosX = repmat(axWH(1)*(0:cols-1),rows,1);
    h4leg = zeros(xCt,1);
    
    ymaxmin = zeros(xCt,2);
    xmaxmin = zeros(xCt,1);
    for iterE = 1:xCt
        if ~dataRef(iterE)
            continue
        end
        if leaveSpaces4omits == 1
            data = grpTable.returnData{sum(dataRef(1:iterE))};
        else
            data = grpTable.returnData{:};
        end
        if isempty(data)
            continue
        end
        ymaxmin(iterE,:) = [max(data(:)) min(data(:))];
        xmaxmin(iterE) = size(data,2);
    end
    xtickpos = linspace(1,max(xmaxmin),5)-1;
    for iterE = 1:xCt
        clearvars plotCt
        if leaveSpaces4omits == 1
            ndxName = grpIDlist{sum(dataRef(1:iterE))};
        else
            ndxName = grpIDlist{iterE};
        end
        ms2frm = @(x) x*(grpTable.recordRate(ndxName)/1000);
        frm2ms = @(x) x/(grpTable.recordRate(ndxName)/1000);
        if explicitX(1) == 1
            minX = ms2frm(explicitX(2));
            maxX = ms2frm(explicitX(3));
            xtickpos = linspace(minX,maxX,explicitX(4));
        end
        
        data = grpTable.returnData{ndxName};
        ytickval = linspace(0,100,5);
        if ~isempty(strfind(graphName,'_TO'))
            data = data(:,4);
            perData = sum(grpTable.jumpTest{ndxName});
        elseif ~isempty(strfind(graphName,'_MVNT'))
            data = data(:,5);
            perData = sum(grpTable.moveTest{ndxName});
        elseif binaryMovement == 0
            if size(data,2) < range(xtickpos)+1 && padDir == 1
                data = cat(2,data,NaN(size(data,1),range(xtickpos)-size(data,2)-1));
            elseif size(data,2) < range(xtickpos)+1 && padDir == -1
                data = cat(2,NaN(size(data,1),range(xtickpos)-size(data,2)-1),data);
            elseif size(data,2) > range(xtickpos)+1 && padDir == 1
                data(:,range(xtickpos)+1:end) = [];
            elseif size(data,2) > range(xtickpos)+1 && padDir == -1
                data(:,1:size(data,2)-range(xtickpos)-1) = [];
            end
            if ~isempty(strfind(graphName,'Accel'))
                ytickval = linspace(-150,150,4);
            elseif ~isempty(strfind(graphName,'Vel'))
                if ~isempty(strfind(graphName,'turning'))
                    if ~isempty(strfind(graphName,'absMotion'))
                        ytickval = linspace(0,360,7);
                    elseif medianAndBounds == 1
                        ytickval = linspace(-0.05,0.05,11);
                    else
                        ytickval = linspace(-0.2,0.2,5)*1;
                    end
                elseif ~isempty(strfind(graphName,'absMotion')) || ~isempty(strfind(graphName,'sidepass')) || ~isempty(strfind(graphName,'movement'))
                    if strcmp(grpTable.stimType{ndxName},'Photoactivation')
                        ytickval = linspace(0,0.5,3);
                    elseif ~isempty(strfind(graphName,'absMotion')) || ~isempty(strfind(graphName,'movement'))
                        ytickval = linspace(0,0.05,3);
                    else
                        ytickval = linspace(-0.5,0.5,3);
                    end
                elseif strcmp(grpTable.stimType{ndxName},'Photoactivation')
                    ytickval = linspace(-0.2,0.2,3);
                else
                    ytickval = linspace(-0.1,0.1,3);
                end
            elseif ~isempty(strfind(graphName,'Dist'))
                ytickval = linspace(-0.5,0.5,3);
                if ~isempty(strfind(graphName,'sidepass')) && strcmp(grpTable.stimType{ndxName},'Photoactivation')
                    ytickval = linspace(0,0.5,3);
                elseif ~isempty(strfind(graphName,'turning'))
                    if strcmp(grpTable.stimType{ndxName},'Photoactivation')
                        ytickval = linspace(0,20,3);
                    else
                        if medianAndBounds == 1
                            ytickval = linspace(-90,90,7);
                        else
                            ytickval = linspace(-90,90,7);
                        end
                        if ~isempty(strfind(graphName,'absMotion'))
                            data = abs(data);
                            ytickval = linspace(0,90,4);
                        end
                    end
                elseif ~isempty(strfind(graphName,'movement')) || ~isempty(strfind(graphName,'absMotion'))
                    ytickval = linspace(0,1,3);
                end
            end
        end
        if explicitY(1) == 1
            ytickval = linspace(explicitY(2),explicitY(3),explicitY(4));
        end
        xAxisName = 'ms';
        ytickpos = linspace(0,1,numel(ytickval));
        if ~dataRef(iterE)
            if iterE == 1
                set(hax(1),'xcolor','none','ycolor','none')
            end
            continue
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
            set(hax(iterE),'xcolor','none','ycolor','none')
            continue
        end
        
        xbaseX = [xtickpos;xtickpos;xtickpos];
        xbaseY = [zeros(size(xtickpos))-0.02;zeros(size(xtickpos))-0.05;zeros(size(xtickpos))-0.02];
        plot(xbaseX(:),xbaseY(:),'linewidth',1,'color',lineCb)
        gridLinesY = repmat(ytickpos(2:end),3,1);
        gridLinesY(3,:) = NaN;
        gridLinesX = repmat([xtickpos(1);xtickpos(end);NaN],1,size(gridLinesY,2));
        plot(gridLinesX(:),gridLinesY(:),'linewidth',1,'color',lineCb)
        if showCountTrace == 1
            countY = sum(~isnan(data),1);
            countYplot = countY/max(countY);
            countX = (xtickpos(1):xtickpos(end));
            plot(countX,countYplot,'color',[0 0 0]+0.5,'linestyle','--','linewidth',1)
        end
        if smoothData == 1
            frms_per_ms = round(grpTable.recordRate(ndxName)/1000);
            ms2frm = @(x) x*frms_per_ms;
            smthLeg = ms2frm(5);
            smthWin = smthLeg*2+1;
            parfor iterS = 1:size(data,1)
                data2smooth = data(iterS,:);
                numRefs = find(~isnan(data2smooth));
                if isempty(numRefs), continue, end
                if numel(numRefs) < smthWin
                    data(iterS,:) = NaN;
                    continue
                end
                data2smoothNums = data2smooth(numRefs);
                data2smoothNums = smooth(data2smoothNums,smthWin);
                data2smoothNums(1:smthLeg) = NaN;
                data2smoothNums(end-smthLeg:end) = NaN;
                data2smooth(numRefs) = data2smoothNums;
                data(iterS,:) = data2smooth(:)';
            end
        end
        if medianAndBounds == 1
            if errorMode(2) >= 2
                segCt = floor(size(data,2)/errorBinSize);
                data = data(:,1:segCt*errorBinSize);
                origSize = size(data);
                data = reshape(data,[size(data,1) errorBinSize segCt]);
                if errorMode(2) == 3
                    dataMax = zeros(origSize(1),1,size(data,3));
                    for iterM = 1:origSize(1)
                        iterData = data(iterM,:,:);
                        for iterN = 1:segCt
                            dataExtreme = [max(iterData(1,:,iterN),[],2) min(iterData(1,:,iterN),[],2)];
                            [~,dataExtremeRef] = max(abs(dataExtreme));
                            dataExtreme = dataExtreme(dataExtremeRef);
                            dataMax(iterM,1,iterN) = dataExtreme;
                        end
                    end
                else
                    dataMax = mean(data,2);
                end
                data = repmat(dataMax,[1 errorBinSize 1]);
                data = permute(data,[2 3 1]);
                data = reshape(data,fliplr(origSize))';
            end
            if errorMode(1) == 1
                meanVec = nanmean(data);
                errA = nanstd(data);
                data = [meanVec+errA;meanVec;meanVec-errA];
            elseif errorMode(1) == 2
                q2 = nanmedian(data,1);
                q1 = prctile(data,25,1);
                q3 = prctile(data,75,1);
%                 errA = 1.57*(q3-q1)/sqrt(numel(meanVec));
                data = [q3;q2;q1];
            end
        elseif medianAndBounds == 2
            iqrVals = NaN(size(data,1),1);
            for iterD = 1:size(data,1)
                dataVec = data(iterD,:);
                if ~isempty(strfind(graphName,'Vel'))
                    iqrVals(iterD) = nanmean(dataVec);
                else
                    dataVec(isnan(dataVec)) = [];
                    if isempty(dataVec), continue, end
                    iqrVals(iterD) = dataVec(end);
                end
            end
            iqrTest = iqrVals > prctile(iqrVals,25) & iqrVals < prctile(iqrVals,75);
            data = data(iqrTest,:);
            plotCt = size(data,1);
        end
        if singleAx == 0
            color = colorOps{1};
        else
            color = colorOps{iterE};
        end
        xplot = [(xtickpos(1):xtickpos(end)) NaN];
        medY = [];
        if binaryMovement == 1
            data(data >= 0.025) = 1;
            data(data < 0.025) = 0;
            yplot = nansum(data)./sum(~isnan(data))*100;
        elseif ~isempty(strfind(graphName,'sum')) || ~isempty(strfind(graphName,'raster'))
            perData = perData-sum(isnan(data));
            yplot = zeros(1,numel(xplot)-1);
            rasterplot = NaN(1,numel(xplot)-1);
            if xScale == 1
                compareVec = linspace(log10(1),log10(frm2ms(xplot(end))),numel(xplot));
                data = log10(data);
            else
                compareVec = frm2ms(xplot);
            end
            for iterD = 1:numel(yplot)
                yplot(iterD) = sum(data <= compareVec(iterD));
                if iterD > 1 && yplot(iterD) > yplot(iterD-1)
                    rasterplot(iterD) = 1;
                    rasterplot(iterD-1) = 0;
                end
            end
            yplot = yplot/perData;
            medY = find(yplot >= 0.5,1,'first');
            if ~isempty(strfind(graphName,'raster'))
                yplot = rasterplot;
            end
        elseif ~isempty(strfind(graphName,'hist'))
            barwidth = ms2frm(barWidth_ms);
            xhist = (xtickpos(1):barwidth:xtickpos(end));
            yhist = hist(ms2frm(data),xhist);
            yplot = interp1(xhist,yhist,xplot(1:end-1));
            yplot = yplot/range(yplot);
            
            fitdata = ms2frm(data);
            fitdata(isnan(fitdata)) = [];
            orig_state = warning;
            warning('off','all')
            xmin = xtickpos(1);      xmax = xtickpos(end);
            fitX = linspace(xmin,xmax,1000);
            gaussFitCt = curveFit(2);
            AIC = Inf(1,gaussFitCt);
            BIC = Inf(1,gaussFitCt);
            ComponentProportion = cell(1,gaussFitCt);
            GMModels = cell(1,gaussFitCt);
            options = statset('MaxIter',500);
            for k = 1:gaussFitCt
                try
                GMModels{k} = fitgmdist(fitdata,k,'Options',options,'CovarianceType','diagonal');
                AIC(k)= GMModels{k}.AIC;
                BIC(k)= GMModels{k}.BIC;
                ComponentProportion(k) = GMModels{k}.ComponentProportion;
                catch
                end
            end
            [~,numComponentsB] = min(BIC);
            [~,numComponentsA] = min(AIC);
            numComponents = ceil(mean([numComponentsA numComponentsB]));
            obj = GMModels{numComponents};
            if ~isempty(obj) && obj.Converged
                while range(obj.ComponentProportion) > 0.75 && numComponents > 1
                    numComponents = numComponents-1;
                    obj = GMModels{numComponents};
                end
                while obj.NegativeLogLikelihood < 200 && numComponents > 1
                    numComponents = numComponents-1;
                    obj = GMModels{numComponents};
                end
            end
            if curveFit(1) == 0
                obj = [];
            end
            warning(orig_state)
        else
            yplot = data;
            yplot = (yplot-ytickval(1))/range(ytickval);
        end
        yplot = [yplot NaN(size(yplot,1),1)]';
        xplot = repmat(xplot(:),1,size(yplot,2));
        jumpVnot = 0;
        if medianAndBounds > 0
            if showCountTrace == 1
                yplot(countY <= 5,:) = NaN;
            end
            plot([xplot(1:end-1,1);flipud(xplot(1:end-1,3));xplot(1,1)],...
                [yplot(1:end-1,1);flipud(yplot(1:end-1,3));yplot(1,1)],...
                'LineWidth',0.5,'color',color);
            hPlot = plot(xplot(:,2),yplot(:,2),'LineWidth',2,'color',color);
        elseif jumpVnot == 0
            if size(yplot,2) > 100
                randNdcs = randperm(size(data,1));
                yplot = yplot(:,randNdcs);
                xplot = xplot(:,randNdcs);
            end
            if ~isempty(strfind(graphName,'hist'))
                yhist = (yhist-min(yhist))/range(yhist);
                hPlot = bar(xhist,yhist,'stack','barwidth',1,'LineWidth',0.5,'facecolor',color.^0.75);
                if ~isempty(obj) && obj.Converged
                    fitY = pdf(obj,fitX');
                    fitY = ((fitY-min(fitY(:)))/range(fitY));
                    plot(fitX,fitY,'color','r','LineWidth',1);
                end
                plot(ms2frm([-75 -75]),[0 1],':','color','r','linewidth',2)
            else
                hPlot = plot(xplot,yplot,'LineWidth',0.5,'color',color);
            end
        else
            for iterBrk = 1:size(yplot,2)
                jumpC = (grpTable.jumpTest{ndxName}(iterBrk)+3);
                hPlot = plot(xplot(:,iterBrk),yplot(:,iterBrk),'LineWidth',1,'linestyle','-','color',colorOps{jumpC});
                hPlot.Color(4) = 0.5;
            end
        end
        if ~isempty(medY) && showMedian == 1
            plot([xplot(medY) xplot(medY)],[0 1],'color','k',...
                'linestyle','--','linewidth',2);
        end
        if singleAx == 1
            h4leg(iterE) = hPlot(1);
        end
        ycolor = fontC;
        xtickval = round(xtickpos/grpTable.recordRate(ndxName)*1000*10)/10;
        if abs(ytickval(1)) == abs(ytickval(end))
            plot([xtickpos(1) xtickpos(end)],[0.5 0.5],'linewidth',1,'color',fontC)
        end
        
        
        if iterE == 1 || ~singleAx
            haxPos = get(hax(iterE),'position');
            haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
                haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
            haxPos = [haxPos(1:2)+[startPosX(iterE) startPosY(iterE)] axWH.*[0.85 0.75]];
            
            set(hax(iterE),'position',haxPos,'color','none','box','off')
            %     continue
            set(hax(iterE),'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1) xtickpos(end)],...
                'YLim',[ytickpos(1)-range(ytickpos)*0.1 ytickpos(end)],'tickdir','out','xtick',[],'yticklabel',[],'ytick',ytickpos,'ycolor','none','xcolor','none')
            for iterL = 1:numel(ytickval)
                text(xtickpos(1)-range(xtickpos)*0.05,ytickpos(iterL),num2str(ytickval(iterL)),'rotation',0,...
                    'color',ycolor,'horizontalalignment','right','fontsize',fontSb);
            end
        end
        
        xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.2];
        if singleAx == 1
            xlabelPos(2) = xlabelPos(2)-(iterE)*range(ytickpos)*0.2;
        end
        if isempty(strfind(graphName,'departureScatter'))
            htxtA = text(xtickpos(1),xlabelPos(2),strtrim(grpTable.plotID{ndxName}),...
                'horizontalalignment','left','rotation',0,...
                'color',fontC,'interpreter','none','fontsize',fontSc);
            jumpCt = sum(grpTable.jumpTest{ndxName});
            if ~exist('plotCt','var')
                plotCt = grpTable.plotCount(ndxName);
            end
            if strcmp(subset2plot,'jumping')
                dataCountTotal = ['n = ' num2str(jumpCt)];
            elseif strcmp(subset2plot,'nonjumping')
                dataCountTotal = ['n = ' num2str(plotCt)];
            else
                dataCountTotal = ['n = ' num2str(plotCt),'    ' num2str(jumpCt),...
                    ' jumping    ' num2str(plotCt-jumpCt) ' nonjumping'];
            end
        else
            dataCountTotal = ['n = ' num2str(sum(grpTable.jumpCount(:)))];
            ylabelPos = [xtickpos(1)-range(xtickpos)*0.2 mean(ytickpos) 1];
            htxtA = text(ylabelPos(1),ylabelPos(2),'Departure Azimuth (degrees)',...
                'horizontalalignment','center','rotation',90,...
                'color',fontC,'interpreter','none','fontsize',fontSa);
        end
        xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.25];
        if singleAx == 1
            xlabelPos(2) = xlabelPos(2)-(iterE)*range(ytickpos)*0.2;
        end
        htxtB = text(xtickpos(1),xlabelPos(2),dataCountTotal,...
            'horizontalalignment','left','rotation',0,'VerticalAlignment','top',...
            'color',fontC,'interpreter','none','fontsize',fontSb,'parent',hax(iterE));
        
        if iterE == xCt
            if ~singleAx
                htxtA.Position(2) = htxtA.Position(2)-range(ytickpos)*0.25;
                htxtB.Position(2) = htxtB.Position(2)-range(ytickpos)*0.25;
            end
            set(get(hax(iterE),'xlabel'),'units','normalized')
            xlabelPos = get(get(hax(iterE),'xlabel'),'position');
            xlabelPos(1:2) = [mean(xtickpos) ytickpos(1)-range(ytickpos)*0.2];
            if xScale == 1
                xtickval = 10.^(linspace(0,log10(xtickval(end)),numel(xtickval)));
            end
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
    if singleAx == 1
        hLeg = legend(h4leg,grpIDlist,'interpreter','none');
        legPos = get(hLeg,'pos');
        legPos(1) = legPos(1)*1.3;
        set(hLeg,'pos',legPos)
    end
    text(xtickpos(1),range(ytickpos)*0.5+ytickpos(end),pdfName,'parent',hax(1),...
        'horizontalalignment','left','interpreter','none',...
        'rotation',0,'color',fontC,'fontsize',round(fontSa*1));
    exportVisualizationFigure(pdfName,saveFig)
end