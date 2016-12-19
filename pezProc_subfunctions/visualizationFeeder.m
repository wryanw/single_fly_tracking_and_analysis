function visualizationFeeder
% clear all

%%%%% computer and directory variables and information
op_sys = system_dependent('getos');
if strfind(op_sys,'Microsoft Windows 7')
    archDir = [filesep filesep 'tier2' filesep 'card'];
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
analysisDir = fullfile(archDir,'Data_pez3000_analyzed');

[~,localUserName] = dos('echo %USERNAME%');
localUserName = localUserName(1:end-1);
repositoryName = 'pezAnalysisRepository';
repositoryDir = fullfile('C:','Users',localUserName,'Documents',repositoryName);
subfun_dir = fullfile(repositoryDir,'pezProc_subfunctions');
saved_var_dir = fullfile(repositoryDir,'pezProc_saved_variables');
addpath(subfun_dir,saved_var_dir)
addpath(fullfile(repositoryDir,'Pez3000_Gui_folder','Matlab_functions','Support_Programs'))

exptIDlist = dir(analysisDir);
exptIDexist = cell2mat({exptIDlist(:).isdir});
exptIDlengthTest = cellfun(@(x) numel(x) == 16,{exptIDlist(:).name});
exptIDlist = {exptIDlist(min(exptIDexist,exptIDlengthTest)).name};

%%

% %%%%% run a single collection
% keepTest = cellfun(@(x) strcmp('0024',x(1:4)),exptIDlist) | cellfun(@(x) strcmp('0047',x(1:4)),exptIDlist);
% exptIDlist = exptIDlist(keepTest)';
% saveDir = 'Y:\WryanW\hack4gwynethrenewal';
% saveName = 'graphingVars_DN_activationAndSilencing';
% runMode = 2;

%%%%% OL Screen
% sheetName = 'exptIDlist';
% excelDir = '\\tier2\card\WryanW\flyPez_experiment_management';
% excelName = 'OL_screen.xlsx';

% exptIDtable = cell2table(exptIDlist);
% exptIDtable.Properties.VariableNames = {'exptID'};
% writetable(exptIDtable,fullfile(excelDir,excelName),'Sheet',sheetName)

% %%%% ming paper
% keepTest = cellfun(@(x) strcmp('0058',x(1:4)),exptIDlist);
% exptIDlist = exptIDlist(keepTest)';
% exptIDlist = [exptIDlist;{'0019000004280129';'0052000003440377';'0052000003490377'
%     '0052000003210377';'0052000020390377';'0052000004300377';'0052000020390369'}];
% sheetName = 'exptIDlist';
% excelDir = 'Y:\Ming_RubinLab';
% excelName = 'mingPaperFinalList.xlsx';
% writetable(cell2table(exptIDlist),fullfile(excelDir,excelName),'Sheet',sheetName)
%%
% excelDir = 'Y:\Ming_RubinLab';
% excelName = 'Table_LC_v2-1.xlsx';
% sheetName = 'Sheet1';
%%
%%%% all DL
exptIDlist = exptIDlist(cellfun(@(x) strcmp('00000430',x(5:12)),exptIDlist))';
excelDir = '\\tier2\card\WryanW\flyPez_experiment_management';
excelName = 'DL_psychophys.xlsx';
sheetName = 'change_this_name';

%%%% methods paper
% keepTestA = cellfun(@(x) strcmp('00000430',x(5:12)),exptIDlist);
% keepTestB = cellfun(@(x) strcmp('0048',x(1:4)),exptIDlist);
% keepTestC = cellfun(@(x) strcmp('0052000004300377',x),exptIDlist);
% keepTestD = cellfun(@(x) strcmp('0377',x(13:end)),exptIDlist);
% exptIDlist = exptIDlist((keepTestA | keepTestB | keepTestD) & ~keepTestC);
% exptIDlist = exptIDlist(:);
% exptIDlist = [exptIDlist;{'0019000003210129';'0019000003440129';'0019000003490129'}];
% exptIDtable = cell2table(exptIDlist);
% exptIDtable.Properties.VariableNames = {'exptID'};
% writetable(exptIDtable,fullfile(excelDir,excelName),'Sheet',sheetName)
% sheetName = 'exptIDlist';
% excelDir = '\\tier2\card\WryanW\flyPez_experiment_management';
% excelName = 'bulk_testing_geno.xlsx';
% 
% excelTable = readtable(fullfile(excelDir,excelName),'Sheet',sheetName);
% exptIDlist = excelTable.exptID;

% %% %%% ryan hack to only run select expt id's from excel file
% filedir = 'Y:\WryanW\flyPez_experiment_management';
% filename = 'bulk_testing_geno.xlsx';
% filepath = fullfile(filedir,filename);
% excelData = dataset('XLSFile',filepath,'Sheet','figure_summary_group');
% includeVec = logical(excelData.include);
% colVec = excelData.Collection_ID(includeVec);
% genoVec = excelData.Genotype_ID(includeVec);
% protoVec = excelData.Protocol_ID(includeVec);
% bulkCt = min([numel(colVec) numel(genoVec) numel(protoVec)]);
% exptIDlist = cell(bulkCt,1);
% for iterB = 1:bulkCt
%     colStr = sprintf('%04s',num2str(colVec(iterB)));
%     genoStr = sprintf('%08s',num2str(genoVec(iterB)));
%     protoStr = sprintf('%04s',num2str(protoVec(iterB)));
%     exptIDlist{iterB} = [colStr,genoStr,protoStr];
% end
% saveDir = 'Z:\CARD LAB SHARE\Ming_temp';
% saveName = 'LC_activation_with_controls';

%%%%%

% guiVarDir = fullfile(dm11Dir,'Pez3000_Gui_folder','Gui_saved_variables');
% groupData = load(fullfile(guiVarDir,'Saved_Group_IDs.mat'));
% groupData = groupData.Saved_Group_IDs;
% [groupNames,refsFirst,refsFull] = unique(groupData.Group_Desc);
% groupNames = strtrim(groupNames)
% groupUsers = groupData.User_ID(refsFirst);
% groupExpts = cell(numel(refsFirst),1);
% for iterGrp = 1:numel(refsFirst)
%     groupExpts(iterGrp) = {groupData.Experiment_ID(refsFull == iterGrp)};
% end
% groupData = [groupNames,groupUsers,groupExpts];
% groupTableData = groupData(:,1:2);
% % exptIDlist = groupData(strcmp(groupTableData(:,1),'CCX_activation_group'),3);
% exptIDlist = groupData(strcmp(groupTableData(:,1),'OL_lines_collection_19_chrimson_only'),3);
% exptIDlist = exptIDlist{1}

%%%% single id
% exptIDlist = {'0044000004300395'};
%%
showDiscards = 0;
infoseek = repmat({''},numel(exptIDlist),1);
typeOps = {'loom','activation'};
typeChoice = 1;
type = typeOps{typeChoice};
for iterA = 1:numel(exptIDlist)
    if iterA == 1
        if strcmp(type,'loom')
            varNames = {'elevation','azimuth','startSize','stopSize','lv','contrast','dataCt'};
        else
            varNames = {'stiminfoA','stiminfoB','stiminfoC','dataCt'};
        end
        addlVars = {'Record_Rate'
            'User_ID'
            'Food_Type'
            'ParentA_name'
            'ParentA_genotype'
            'ParentB_name'
            'ParentB_genotype'
            'ParentA_ID'
            'ParentB_ID'
            'Males'
            'Females'}';
        varNames = cat(2,varNames,{'Synonyms'},addlVars);
        exptIDtable = cell2table(cell(numel(exptIDlist),numel(varNames)),'RowNames',exptIDlist,...
            'VariableNames',varNames);
    end
    exptID = exptIDlist{iterA};
    
        expt_results_dir = fullfile(analysisDir,exptID);
        exptInfoMergedName = [exptID '_experimentInfoMerged.mat'];
        exptInfoMergedPath = fullfile(expt_results_dir,exptInfoMergedName);
        if ~exist(exptInfoMergedPath,'file')
            exptIDtable(exptID,:) = [];
            if showDiscards == 1
                disp([exptID ' - missing expt info'])
            end
            continue
        end
        experimentInfoMerged = load(exptInfoMergedPath,'experimentInfoMerged');
        exptInfo = experimentInfoMerged.experimentInfoMerged;

        graphTablePath = fullfile(expt_results_dir,[exptID '_dataForVisualization.mat']);
        if ~exist(graphTablePath,'file')
            exptIDtable(exptID,:) = [];
            if showDiscards == 1
                disp([exptID ' - no data for visualization file'])
            end
            continue
        end
        graphTableLoading = load(graphTablePath);
        graphTable = graphTableLoading.graphTable;
        dataCt = numel(strcmp(graphTable.finalStatus,'Graphable'));
        exptIDtable.dataCt(exptID) = {dataCt};
        
    infoseek{iterA} = exptInfo.Stimuli_Type{1};
    stimInfo = exptInfo.Stimuli_Type{1};
    splitCell = strsplit(stimInfo,'_');
    if strcmp(type,'loom')
        if ~strcmp(splitCell{1},'loom')
            exptIDtable(exptID,:) = [];
            if showDiscards == 1
                disp([exptID ' - expt type mismatch'])
            end
            continue
        end
        exptIDtable.elevation(exptID) = {exptInfo.Stimuli_Vars.Elevation};
        exptIDtable.azimuth(exptID) = {exptInfo.Stimuli_Vars.Azimuth};
        exptIDtable.lv(exptID) = {splitCell{3}(strfind(splitCell{3},'lv')+2:end)};
        exptIDtable.startSize(exptID) = {splitCell{2}(1:strfind(splitCell{2},'to')-1)};
        exptIDtable.stopSize(exptID) = {splitCell{2}(strfind(splitCell{2},'to')+2:end)};
        exptIDtable.contrast(exptID) = splitCell(end);
    else
        if strcmp(splitCell{1},'loom')
            exptIDtable(exptID,:) = [];
            if showDiscards == 1
                disp([exptID ' - expt type mismatch'])
            end
            continue
        end
        if ~iscell(exptInfo.Photo_Activation{1})
            exptInfo.Photo_Activation = {exptInfo.Photo_Activation};
        end
%         testA = strcmp(exptInfo.Photo_Activation{1}(1),'pulse_General_widthBegin50_widthEnd50_cycles1_intensity50');
%         testB = strcmp(exptInfo.Photo_Activation{1}(1),'pulse_General_widthBegin50_widthEnd50_cycles2_intensity50');
%         if ~(testA || testB)
%             exptIDtable(exptID,:) = [];
%             continue
%         end
        for iterS = 1:numel(exptInfo.Photo_Activation{1})
            exptIDtable.(varNames{iterS+2})(exptID) = exptInfo.Photo_Activation{1}(iterS);
        end
    end
    for iterAddl = 1:numel(addlVars)
        addlInfo = exptInfo.(addlVars{iterAddl});
        while iscell(addlInfo)
            addlInfo = addlInfo{1};
        end
        exptIDtable.(addlVars{iterAddl})(exptID) = {addlInfo};
    end
end
%% identify synonymous experiments and save
fields2ignore = {'Synonyms','dataCt','User_ID'};
exptIDtable.Synonyms = cell(size(exptIDtable,1),1);
exptIDlist = exptIDtable.Properties.RowNames;
for iterA = 1:numel(exptIDlist)
    nextSyn = max(cell2mat(exptIDtable.Synonyms));
    if isempty(nextSyn)
        nextSyn = 1;
    else
        nextSyn = nextSyn+1;
    end
    exptIDtable.Synonyms(iterA) = {nextSyn};
    for iterB = 1:numel(exptIDlist)
        if iterB == iterA, continue, end
        testTableA = exptIDtable(iterA,:);
        testTableB = exptIDtable(iterB,:);
        for iterI = 1:numel(fields2ignore)
            testTableA.(fields2ignore{iterI}) = [];
            testTableB.(fields2ignore{iterI}) = [];
        end
        testTableA = table2cell(testTableA);
        testTableB = table2cell(testTableB);
        if isequal(testTableA,testTableB)
            if ~isempty(exptIDtable.Synonyms{iterB})
                exptIDtable.Synonyms(iterA) = exptIDtable.Synonyms(iterB);
            end
        end
    end
end
rowNames = exptIDtable.Properties.RowNames;
varNames = exptIDtable.Properties.VariableNames;
writeIDtable = table2cell(exptIDtable);
rowNames = cellfun(@(x) ['ID#' x],rowNames,'uniformoutput',false);
writeIDtable = cat(2,rowNames,writeIDtable);
writeIDtable = cell2table(writeIDtable,'VariableNames',[{'exptID'} varNames]);
writetable(writeIDtable,fullfile(excelDir,excelName),'Sheet',sheetName,'WriteRowNames',true)
%% labels for looming experiments
readTable = readtable(fullfile(excelDir,excelName),'Sheet','exptIDtable','ReadRowNames',true);
readTable.labelA = cell(size(readTable,1),1);
rowNames = readTable.Properties.RowNames;
querryName = {'azimuth','elevation','startSize','stopSize','lv','contrast'};
shortName = {'azi','ele','strt','stop','lv',''};
for iterL = 1:numel(rowNames)
    querryInfo = cellfun(@(x) readTable.(x)(rowNames{iterL}),querryName(1:end-1),'uniformoutput',false);
    querryInfo = strtrim(cellstr(num2str(cell2mat(querryInfo'))))';
    querryInfo = cat(2,querryInfo,readTable.(querryName{end})(rowNames{iterL}));
    querryInfo = cat(1,shortName,querryInfo,repmat({'_'},1,numel(querryInfo)));
    querryInfo = querryInfo(:)';
    querryInfo = cat(2,querryInfo{1:end-1});
    querryInfo = regexprep(querryInfo,'_blackonwhite','');
    querryInfo = regexprep(querryInfo,'_whiteonwhite','_whiteStim');
    readTable.labelA(rowNames{iterL}) = {querryInfo};
end
if numel(unique(readTable.Synonyms)) ~= numel(unique(readTable.labelA))
    disp('unique label count does not match unique experiment count')
end
writetable(readTable,fullfile(excelDir,excelName),'Sheet','exptIDtable','WriteRowNames',true)
