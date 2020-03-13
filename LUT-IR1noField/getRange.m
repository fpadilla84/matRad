function r = getRange(z, dose, doPlot)
if nargin < 3
    doPlot = false;
end
% doPlot = true;
[maxDose, indexMaxDose] = max(dose);
normalizedToMaxDdose = dose./maxDose;

[~, index50percent] = getPositionOfDoseLevel(z, dose, 0.5);

zPeakTo50Percent = z(indexMaxDose:index50percent);
dPeakTo50Percent = normalizedToMaxDdose(indexMaxDose:index50percent);
if numel(zPeakTo50Percent) < 3
    degreeFit = 1;
else
    degreeFit = 2;
end
ws = warning('off','all');  % Turn off warning

p2 = polyfit(zPeakTo50Percent(:), dPeakTo50Percent(:), degreeFit);
warning(ws);  % Turn it back on.
zTempPeakTo50Percent = zPeakTo50Percent(1):0.001:zPeakTo50Percent(end);
[r,~] = getPositionOfDoseLevel(zTempPeakTo50Percent, polyval(p2, zTempPeakTo50Percent), 0.8);
if doPlot
    figure
    plot(zPeakTo50Percent, dPeakTo50Percent,'.', zTempPeakTo50Percent, polyval(p2, zTempPeakTo50Percent))
end
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

