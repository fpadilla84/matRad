function [cumWeightFinalFromFile,cumulativeMetersetWeight ]= readGateTPSsourceF(theDir, fileName, mIsFastLookUp)
% myDict=containers.Map;

% theDir='/home/aresch/simulations/RS_LETcalc/adrian_9thJuly2017/outTemp/';
% fileName='SOBP_MC_deg100.txt';
if nargin <3
    mIsFastLookUp=true;
end
fid=fopen(fullfile(theDir,fileName));
k=0;
tline = fgetl(fid);
beginWeights=false;beginNumberOfControlPoints=false;
beginCumWeightFinalFromFile=false;
sumControlPoints=1;
cumulativeMetersetWeight=0;controlPointsFile=-1;cumWeightFinalFromFile=-1;
while ischar(tline)
    k=k+1;
%     disp(tline)
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    
    %% cum Weights
    if beginCumWeightFinalFromFile
        numCon=regexp(tline,'([\-\+0-9.eE]+)','tokens');
        cumWeightFinalFromFile=str2num(numCon{1}{end});
        beginCumWeightFinalFromFile=false;
    end
    if ~isempty(regexp(tline,'[# \t]+(TotalMetersetWeightOfAllFields)','Once'))
         beginCumWeightFinalFromFile=true;
    end
    if ~mIsFastLookUp
        %% weights
        if beginWeights
            xyw_line=regexp(tline,'([\-\+0-9.]+)[ \t]+([\-\+0-9.]+)[ \t]+([\-\+0-9.e]+)','tokens');
            if ~isempty(xyw_line)
                cumulativeMetersetWeight=cumulativeMetersetWeight+str2double(xyw_line{1}{3});
            else

                beginWeights=false;
                sumControlPoints=sumControlPoints+1;
                if mod(sumControlPoints,10)==0
                    disp(tline)
                    disp(sumControlPoints)
                end
            end
        end
        if ~isempty(regexp(tline,'[# \t]+[Xx \tYu]+(Weight)','Once'))
            beginWeights=true;
        end
    end
    %% control points
    if ~beginWeights && beginNumberOfControlPoints
        numCon=regexp(tline,'([\-\+0-9.]+)','tokens');
        controlPointsFile=str2num(numCon{1}{end});
        beginNumberOfControlPoints = false;
    end
    if ~beginWeights && ~isempty(regexp(tline,'[# \t]+(NumberOfControlPoints)','Once'))
        beginNumberOfControlPoints=true;
    end
    


%     if k==690
%         break
%     end
end

fclose(fid);

disp(['Control points in file: ', sprintf('%u',controlPointsFile)]);
disp(['Evaluated control points: ', sprintf('%u',sumControlPoints)]);
disp(['Total cumulative meterset weight sum: ', sprintf('%.4e',cumulativeMetersetWeight)]);
disp(['Total cumulative meterset weight from file: ', sprintf('%.4e',cumWeightFinalFromFile)]);
if ~mIsFastLookUp
    disp(['Relative difference: ',sprintf('%.4f',(1-cumWeightFinalFromFile./cumulativeMetersetWeight)*100),'%']);
    if abs(controlPointsFile - sumControlPoints)>0
        error('Control points do not agree: check code and txt file manually')
    end
end
% 
% allItems = regexp(simFileStr,'[# \t]+([a-zA-Z ]+)[ \t]+=[ \t]+([0-9e+-.]+)','tokens');
% if ~isempty(allItems)
%     for i=1:numel(allItems)
%         tempName=regexprep(allItems{i}{1},'[ \t]','');
% %         tempV  =allItems{i}{2};
%         try
%             tempV = str2double(allItems{i}{2});
%         catch
%             tempV = allItems{i}{2};
%         end
%         myDict.(tempName) =tempV;
%     end
% end
% myDict=0;
end