function myDict = readGateStatisticFile(theDir, fileName)
% myDict=containers.Map;
simFileStr=fileread(fullfile(theDir,fileName));
allItems = regexp(simFileStr,'[# \t]+([a-zA-Z ]+)[ \t]+=[ \t]+([0-9e+-.]+)','tokens');
if ~isempty(allItems)
    for i=1:numel(allItems)
        tempName=regexprep(allItems{i}{1},'[ \t]','');
        if ~isempty(tempName)
    %         tempV  =allItems{i}{2};
            try
                tempV = str2double(allItems{i}{2});
            catch
                tempV = allItems{i}{2};
            end
            myDict.(tempName) =tempV;
        end
    end
end
% myDict=0;
end