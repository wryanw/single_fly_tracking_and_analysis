function colorMapVar = colorMapMaker(mapChoice)
%colorMapMaker Makes heatmap using custom inputs, for exploration
%   Detailed explanation goes here
%%
res = 7;
if isempty(mfilename) || ~exist('mapChoice','var')
    cla
    mapChoice = 2;
end
if mapChoice == 0
    colorMapDef = {round([res/4 1 0;NaN 0 0;res/2 0 1])
        round([res/1.7 1 0;NaN 0 0])
        round([NaN 0 0;res/5 0 1;res/5 1 1;res/3 1 0])};
elseif mapChoice == 1
    colorMapDef = {round([res/4 1 0;NaN 0 0;res/2 0 1])
        round([res/1.7 1 0;NaN 0 0])
        round([NaN 0 0;res/5 0 1;res/5 1 1;res/3 1 0])};
elseif mapChoice == 2
    colorMapDef = {[NaN 1 0;NaN 0 0;NaN 0 1]
        [NaN 1 0;NaN 0 0;NaN 0 0]
        [NaN 1 1;NaN 1 0;NaN 0 0]};
else
    colorMapDef = mapChoice;
end
colorMapVar = [];
for i = 1:3
    smallDef = colorMapDef{i};
    if sum(isnan(smallDef(:))) > 0
        smallDef(isnan(smallDef)) = (res-sum(smallDef(~isnan(smallDef(:,1)),1)))/sum(isnan(smallDef(:)));
    end
    channel = [];
    for k = 1:size(smallDef,1)
        legsize = smallDef(k,1);
        legA = smallDef(k,2);
        legB = smallDef(k,3);
        leg = linspace(legA,legB,legsize)';
        channel = cat(1,channel,leg);
    end
    colorMapVar = cat(2,colorMapVar,channel);
end
colormap(colorMapVar)
if isempty(mfilename) || ~exist('colorMapDef','var')
    plot(colorMapVar(:,1),'r')
    hold on
    plot(colorMapVar(:,2),'g')
    plot(colorMapVar(:,3),'b')
    colorbar
end
end

