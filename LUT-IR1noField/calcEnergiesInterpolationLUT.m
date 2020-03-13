% Calculate the required energies for LUT generation 

for i = 1:length(machine.data);
    
    E(i) = machine.data(i).energy .' ;
    peakPos(i) = machine.data(i).peakPos .' ;  
    maxEdep(i) = machine.data(i).maxEdep .' ;
    
end;

[p,ErrorEst] = polyfit(peakPos,E,3) ;
newPos = 20:2:400 ;
newEnergies = polyval(p,newPos,ErrorEst) ;

figure;

plot(peakPos,E,'x',newPos,newEnergies,'-' ) ;

dlmwrite('EnergiesInterpolation.txt', newEnergies .') ;



