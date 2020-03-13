% interpolate all LUT elements to intermediate energies 

% Read calibration data LUT for simulated/measured energies  

availableEnergies = zeros; 
availablePeakPos = zeros; 
availableMaxEdep = zeros; 
energySteps = zeros; 

for i = 1:length(machine.data);
    
    availableEnergies(i) = machine.data(i).energy   ;
        if i== 1. energySteps(i) = machine.data(i+1).energy - machine.data(i).energy ;
        else      energySteps(i) = machine.data(i).energy - machine.data(i-1).energy ;
        end;
    availablePeakPos (i) = machine.data(i).peakPos  ;  
    availableMaxEdep (i) = machine.data(i).maxEdep  ;
    
end;    

% peakPos required for interpolation, energies fixed to values available for irradiation ;

IR1Energies = xlsread('EnergiesIrradiation_IR1MA.xlsx','Energy-RangeDep') ;

newEnergies = IR1Energies(:,1) ;
newPeakPos = interp1(availableEnergies, availablePeakPos, newEnergies) ; 
% newPeakPos = min(availablePeakPos):2:max(availablePeakPos) ;
% [p,ErrorEst] = polyfit(availablePeakPos,availableEnergies,3) ;
% newEnergies = polyval(p,newPeakPos,ErrorEst) ;

% newEnergies = interp1(availablePeakPos,availableEnergies,newPeakPos) ; 
% newEnergies = newEnergies(newEnergies>min(availableEnergies) & newEnergies<max(availableEnergies)) ;
% 
for j=1:length(newEnergies) 
% Energy for new LUT
% Add condition to warranty is between limits of available energies 

newEnergy = newEnergies(j) ;
% find closest energy in simulated LUT array
[Diff index] = min(newEnergy - availableEnergies((newEnergy - availableEnergies)> 0) );

% calculate new IDD function interpolating from closest available LUT  
% scale new Bragg peak position and values 
%     if min(availableEnergies) < newEnergy < max(availableEnergies)
    newMaxPos = interp1(availableEnergies,availablePeakPos,newEnergy) ;
    newMaxEdep = interp1(availableEnergies,availableMaxEdep,newEnergy) ;
    scalingMaxPos = machine.data(index).peakPos ./ newMaxPos ;
    scalingMaxEdep = machine.data(index).maxEdep ./ newMaxEdep ;
%     else 
%     newMaxPos = interp1(availableEnergies,availablePeakPos,max(availableEnergies)) ;
%     newMaxEdep = interp1(availableEnergies,availableMaxEdep,max(availableEnergies)) ;
%     scalingMaxPos = machine.data(index).peakPos ./ newMaxPos ;
%     scalingMaxEdep = machine.data(index).maxEdep ./ newMaxEdep ;    
%     end;
    
% evaluate new IDD function 
    scaledDepths = machine.data(index).depths .* scalingMaxPos ;
    minDepth = min(machine.data(index).depths) ;
    if scaledDepths <= min(machine.data(index).depths)  scaledEdep = interp1(machine.data(index).depths, machine.data(index).Z, min(machine.data(index).depths)) ./ scalingMaxEdep  ; 
    else scaledEdep = interp1(machine.data(index).depths, machine.data(index).Z, scaledDepths) ./ scalingMaxEdep  ;
    end;
    
% Interpolate rest of the parameters from closest available LUT  
scalingFactorEnergy = Diff ./ energySteps(index) ;
newMean1X = machine.data(index).mean1x + scalingFactorEnergy .* (machine.data(index+1).mean1x - machine.data(index).mean1x) ;
newMean2X = machine.data(index).mean2x + scalingFactorEnergy .* (machine.data(index+1).mean2x - machine.data(index).mean2x) ;
newSigma1X = machine.data(index).sigma1x + scalingFactorEnergy .* (machine.data(index+1).sigma1x - machine.data(index).sigma1x) ;
newSigma2X = machine.data(index).sigma2x + scalingFactorEnergy .* (machine.data(index+1).sigma2x - machine.data(index).sigma2x) ;
newfractionAreax = machine.data(index).fractionAreax + scalingFactorEnergy .* (machine.data(index+1).fractionAreax - machine.data(index).fractionAreax) ;

newMean1Y = machine.data(index).mean1y + scalingFactorEnergy .* (machine.data(index+1).mean1y - machine.data(index).mean1y) ;
newMean2Y = machine.data(index).mean2y + scalingFactorEnergy .* (machine.data(index+1).mean2y - machine.data(index).mean2y) ;
newSigma1Y = machine.data(index).sigma1y + scalingFactorEnergy .* (machine.data(index+1).sigma1y - machine.data(index).sigma1y) ;
newSigma2Y = machine.data(index).sigma2y + scalingFactorEnergy .* (machine.data(index+1).sigma2y - machine.data(index).sigma2y) ;
newfractionAreay = machine.data(index).fractionAreay + scalingFactorEnergy .* (machine.data(index+1).fractionAreay - machine.data(index).fractionAreay) ;

newAmpCorrFactor = machine.data(index).AmpCorrFactor + scalingFactorEnergy .* (machine.data(index+1).AmpCorrFactor - machine.data(index).AmpCorrFactor) ;

jdx = j + 15. ;
machine.data(jdx).energy = newEnergy ;
machine.data(jdx).depths = machine.data(index).depths ;
machine.data(jdx).Z = scaledEdep  ;
machine.data(jdx).range = getRange(machine.data(jdx).depths, machine.data(jdx).Z);
machine.data(jdx).peakPos = 0.1 .* getmaxPosition(machine.data(jdx).depths, machine.data(jdx).Z);
machine.data(jdx).maxEdep = max(machine.data(jdx).Z);
machine.data(jdx).BField = machine.data(index).BField ;
machine.data(jdx).offset = machine.data(index).offset ;
machine.data(jdx).initFocus = machine.data(index).initFocus ;
machine.data(jdx).LET = machine.data(index).LET ;

machine.data(jdx).mean1x = newMean1X ;
machine.data(jdx).mean2x = newMean2X ;
machine.data(jdx).sigma1x = newSigma1X ;
machine.data(jdx).sigma2x = newSigma2X ;
machine.data(jdx).fractionAreax = newfractionAreax ;
machine.data(jdx).mean1y = newMean1Y ;
machine.data(jdx).mean2y = newMean2Y ;
machine.data(jdx).sigma1y = newSigma1Y ;
machine.data(jdx).sigma2y = newSigma2Y ;
machine.data(jdx).fractionAreay = newfractionAreay ;
machine.data(jdx).AmpCorrFactor = newAmpCorrFactor ;


end;

% clear all;

    
    