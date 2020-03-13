% close all;
% load machine data
load('protons_IR1noField-GATE-new.mat') ;
% set x,y limits to perform calculation 
Y = -199:2:199 ;
YmidVoxel = -45.5:2:45.5 ;
Z = -199:2:199 ;
ZmidVoxel = -49.5:2:49.5 ;

% figure;
% plot(X,pdfYinX(200,:)); hold on;
% plot(XmidVoxel,pdfYinmidVoxel(200,:)); hold on;


for k = 1:length(machine.data)
k = 7 ; 

LatY = matRad_expGaussFunction(Y,machine.data(k).meanx, machine.data(k).sigmax, machine.data(k).leftGammax, machine.data(k).rightGammax, machine.data(k).fractionleftx, machine.data(k).fractionrightx) ;
LatZ = matRad_expGaussFunction(Z,machine.data(k).meany, machine.data(k).sigmay, machine.data(k).leftGammay, machine.data(k).rightGammay, machine.data(k).fractionlefty, machine.data(k).fractionrighty) ;
%             [xxi,yyi]     = meshgrid(Y,Z);
%             [pdfxx,pdfyy] = meshgrid(LatY,LatZ);
%             % Calculate combined pdf, under assumption of independence
%             pdfxy = pdfxx.*pdfyy; 
% figure;
% % Plot the results
% mesh(xxi,yyi,pdfxy)

sumLatY = sum(LatY,2) ;
sumLatZ = sum(LatZ,2) ;

            for i=1:length(machine.data(1).depths) ;
            % Estimate a continuous pdf from the discrete data
            pdfx = matRad_expGaussFunction(X,machine.data(k).meanx(i), machine.data(k).sigmax(i), machine.data(k).leftGammax(i), machine.data(k).rightGammax(i), machine.data(k).fractionleftx(i), machine.data(k).fractionrightx(i)) ;
            pdfy = matRad_expGaussFunction(Z,machine.data(k).meany(i), machine.data(k).sigmay(i), machine.data(k).leftGammay(i), machine.data(k).rightGammay(i), machine.data(k).fractionlefty(i), machine.data(k).fractionrighty(i)) ;
            % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
%             [xxi,yyi]     = meshgrid(X,Z);
%             [pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
%             % Calculate combined pdf, under assumption of independence
%             pdfxy = pdfxx.*pdfyy; 
            sum2dpdfKernel(i) = sum(sum(pdfxy,2),1) .' ;

            end

% figure;
% Plot the results
% mesh(xxi,yyi,pdfxy)


subplot(7,2,k)
title(strcat('Energy-',num2str(machine.data(k).energy),'mm'))
yyaxis left
plot(10.*machine.data(k).depths, machine.data(k).Z ./ max(machine.data(k).Z(:))) ;
xlim([1 10.*machine.data(k).peakPos + 50]) ;
ylim([0 1.2]) ;

yyaxis right
plot(sumLatY); hold on; plot(sumLatZ) ; hold on;
plot(sum2dpdfKernel) ;
xlim([1 10.*machine.data(k).peakPos + 50]) ;
ylim([0.99 1.01]) ;


end

legend({'LUT-IDD','intLatY', 'intLatZ' , 'int2Dkernel'},'Location','north');

% spot = LatY .* LatZ ;
% dose = machine.data(1).Z * LatY .* LatZ  ;
% 
% 
% figure; 
% subplot(2,2,1)
% contourf(LatY) ; 
% ylim([1 319]) ;
% 
% subplot(2,2,2)
% contourf(LatZ) ; 
% ylim([1 319]) ;
% 
% subplot(2,2,3)
% contourf(spot) ; 
% ylim([1 319]) ;
% 
% subplot(2,2,4)
% contourf(dose) ; 
% ylim([1 319]) ;


% figure; plot(machine.data(1).depths, machine.data(1).sigmax) ;
% figure; plot(machine.data(1).depths, machine.data(1).fractionleftx) ;
% figure; plot(machine.data(1).depths, machine.data(1).fractionrightx) ;


    
    