function r = getRange20(z, dose)
if nargin < 3
    doPlot = false;
end
% doPlot = true;
[maxDose, indexMaxDose] = max(dose);
normalizedToMaxDdose = dose./maxDose;

[r,~] = getPositionOfDoseLevel(z, dose, 0.2);

end

function [ zPos , index ] = getPositionOfDoseLevel( z, D, doseLevel, firstIntersection )
%gives back position and index of a certain dose level. Dose fall off by
%default
%   Detailed explanation goes here
doFirstIntersection = false;
if nargin == 4
    doFirstIntersection = true;
end
D = D./max(D);

D_greater_threshold = D > doseLevel;
indices_D_greater_threshold = find(D_greater_threshold);
if doFirstIntersection
    index = indices_D_greater_threshold(1);
else
    index = indices_D_greater_threshold(end);
end
zPos = z(index);
    


end

