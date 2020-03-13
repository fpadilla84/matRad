Energies = dlmread('dataAvailable.dat');
Energies = round(Energies .* 10) ./ 10 ;


% for i = 1:length(Energies); 
%     
%     currentEnergy = Energies(i); 
%     index = 174. + i ;
%     
%         filename1 = strcat('/home/fpadilla/Matlab15/bin/work/MatRad/LUT-IR1NoField/ParametersYPB-water-', num2str(currentEnergy,'%.1f'), '-IR1noField.dat') ;  % Character Array 
%         filename2 = strcat('/home/fpadilla/Matlab15/bin/work/MatRad/LUT-IR1NoField/ParametersXPB-water-', num2str(currentEnergy,'%.1f'), '-IR1noField.dat') ;  % Character Array 
%         filename3 = strcat('/home/fpadilla/Matlab15/bin/work/MatRad/LUT-IR1NoField/Edep-water-', num2str(currentEnergy,'%.1f'), '-IR1noField.dat') ;  % Character Array 
% 
%     nucLUTy = dlmread(filename1, ' ', 1, 0);
%     nucLUTx = dlmread(filename2, ' ', 1, 0);
%     PDD = dlmread(filename3, ' ', 1, 0);
% 
%     
% machine.data(index).energy = currentEnergy ;
% machine.data(index).BField = 0. ;
% machine.data(index).offset = 0. ;
% 
% machine.data(index).depths = 10 .* PDD(:,1);
% machine.data(index).Z = PDD(:,2);
% 
% machine.data(index).range = 10 .* getRange(PDD(:,1),PDD(:,2));
% machine.data(index).peakPos = getmaxPosition(PDD(:,1),PDD(:,2));   
% 
% 
% machine.data(index).sigma = 10 .* nucLUTx (:,4);
% 
% machine.data(index).meany =  10 .* nucLUTy (:,3);
% machine.data(index).sigmay =  10 .* nucLUTy (:,4);
% machine.data(index).leftGammay = 10 .*  nucLUTy (:,5);
% machine.data(index).rightGammay = 10 .*  nucLUTy (:,6);
% machine.data(index).fractionlefty = nucLUTy (:,7);
% machine.data(index).fractionrighty = nucLUTy (:,8);
% 
% machine.data(index).meanx = 10 .* nucLUTx (:,3);
% machine.data(index).sigmax = 10 .* nucLUTx (:,4);
% 
% machine.data(index).leftGammax = 10 .* nucLUTx (:,5);
% machine.data(index).rightGammax = 10 .* nucLUTx (:,6);
% machine.data(index).fractionleftx = nucLUTx (:,7);
% machine.data(index).fractionrightx = nucLUTx (:,8);
% 
% 
%     
% end;

for i = 1:length(Energies); 
    
    currentEnergy = Energies(i); 
    realEnergies = [62.4 81.3 97.4 111.6 124.7 136.8 148.2 159 169.3 179.2 188.7 198 207 215.7] ;
    index = i+5 ;

        filename1 = strcat('/Users/fpc/Documents/MATLAB/matRad/LUT-IR1noField/Gate-data/doubleGauss-fitting/ParametersZDoubleGauss_noField_', num2str(currentEnergy,'%.0f'), 'MeV.txt') ;  % Character Array 
        filename2 = strcat('/Users/fpc/Documents/MATLAB/matRad/LUT-IR1noField/Gate-data/doubleGauss-fitting/ParametersYDoubleGauss_noField_', num2str(currentEnergy,'%.0f'), 'MeV.txt') ;  % Character Array 
        filename3 = strcat('/Users/fpc/Documents/MATLAB/matRad/LUT-IR1noField/Gate-data/doubleGauss-fitting/dose-PDD-Edep-', num2str(currentEnergy,'%.0f'), 'MeV.txt') ;  % Character Array 

%         filename1 = strcat('/home/fpadilla/MatlabR2017a/matRad/LUT-IR1noField/Gate-data/doubleGauss-fitting/ParametersYDoubleGauss_', num2str(currentEnergy,'%.0f'), 'MeV.txt') ;  % Character Array 
%         filename2 = strcat('/home/fpadilla/MatlabR2017a/matRad/LUT-IR1noField/Gate-data/doubleGauss-fitting/ParametersXDoubleGauss_', num2str(currentEnergy,'%.0f'), 'MeV.txt') ;  % Character Array 
%         filename3 = strcat('/home/fpadilla/MatlabR2017a/matRad/LUT-IR1noField/Gate-data/doubleGauss-fitting/dose-PDD-Edep-', num2str(currentEnergy,'%.0f'), 'MeV.txt') ;  % Character Array 

%         filename1 = strcat('/home/fpadilla/Matlab15/bin/work/MatRad/LUT-IR1BeamModel-Danfysik-fringeFields/GATE-data/nucCorrLatFitX_GATE_B-Danfisyk_IR1BeamModel_', num2str(currentEnergy,'%.0f'), '.dat') ;  % Character Array 
%         filename2 = strcat('/home/fpadilla/Matlab15/bin/work/MatRad/LUT-IR1BeamModel-Danfysik-fringeFields/GATE-data/nucCorrLatFitY_GATE_B-Danfisyk_IR1BeamModel_', num2str(currentEnergy,'%.0f'), '.dat') ;  % Character Array 

    nucLUTy = dlmread(filename1, '\t', 0, 0);
    nucLUTx = dlmread(filename2, '\t', 0, 0);
    scoringDepthsNucLUT = nucLUTy(:,1) ; 
    PDD = dlmread(filename3, ' ', 6, 0);
    PDD = PDD(:,2) .* 100 ./ 10000000 ;
    PDD = flipud(PDD);
%     PDD = dlmread(filename3, ' ');

machine.data(index).energy = realEnergies(i) ;
machine.data(index).BField = 0. ;
machine.data(index).offset = 0. ;

% get PDD with 0.1mm steps
penDepths = 0:1:length(PDD)-1;
penDepths = penDepths .* 0.1  ;
machine.data(index).depths = penDepths .' ;
machine.data(index).range = getRange(penDepths,PDD) ;
machine.data(index).peakPos = getmaxPosition(penDepths,PDD) .* 0.1;  
machine.data(index).maxEdep = max(PDD) ;
machine.data(index).Z = PDD ;


% calculate nucLUTs with 0.1 mm steps
machine.data(index).sigma1 = (interp1(scoringDepthsNucLUT,nucLUTy (:,5),penDepths)) .';
machine.data(index).sigma1(isnan(machine.data(index).sigma1y)) = nucLUTy(1,5) ; 
machine.data(index).sigma2 = (interp1(scoringDepthsNucLUT,nucLUTy (:,8),penDepths)) .';
machine.data(index).sigma2(isnan(machine.data(index).sigma2y)) = nucLUTy(1,8) ; 
machine.data(index).weight = (interp1(scoringDepthsNucLUT, nucLUTy (:,9), penDepths)) .';
machine.data(index).weight(isnan(machine.data(index).fractionAreay)) = nucLUTy (1,9)  ; 


machine.data(index).initFocus = machine.data(1).initFocus ;

    
end;



