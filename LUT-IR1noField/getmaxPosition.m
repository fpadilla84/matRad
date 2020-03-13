function r = getmaxPosition(z, dose)

[maxDose, indexMaxDose] = max(dose);
normalizedToMaxDdose = dose./maxDose;

r = indexMaxDose ;

end
