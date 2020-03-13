function d = gateTPS_doseReScale(edep,theDirPath,gateRTplanFileName)
if ~strcmpi(theDirPath(end),filesep)
    theDirPath=[theDirPath,filesep];
end


    t=readGateStatisticFile([theDirPath,'output'],'Statistics-autoMerged.txt');
    cumulativeMeterSetWeight=readGateTPSsourceF([theDirPath,'data'], gateRTplanFileName,false);
%     doseSTDFromFile = theFMhdSTD.data*sqrt(t.Combinedstatisticsfiles);
    d=edep*(cumulativeMeterSetWeight/t.NumberOfEvents);


end