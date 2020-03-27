close all
clear all

[~,userName]=system('whoami');
if ~isempty(regexp(userName,'fpadilla','Once'))
    tPathGATE='/home/fpadilla/Simulations/IR1-BeamModel/BioSetup/81.3MeV';

else ~isempty(regexp(userName,'fpc','Once'))
     tPathGATE='/Users/fpc/Documents/GATE-Simulations/BioSetup/100.4MeV';
    
end

% Check to make sure that folders actually exists.  Warn user if it doesn't.
if ~isdir(tPathGATE)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', tPathGATE);
  uiwait(warndlg(errorMessage));
  return;
end


%% read GATE 3D dose distributions
mhdF=read_mhd([tPathGATE,'/output'],'phantom_dose3d-Edep-autoMerged.mhd');
% Matrix covering 300 mm x 120 mm (full phantom geometry)
RBE = 1.1 ;
voxelSize = 1 ;
MaterialDensity = 1. ; % Water

dataGATE = mhdF.data ;

dimGATE = size(dataGATE) ;
centraldoseGATE = squeeze(dataGATE(:,round(dimGATE(2)/2),:));
centralSOBPGATE = centraldoseGATE(:,round(dimGATE(3)/2));
numberVoxels = 2 ;

% get mean value of laterally integrated dose YZ in 1 cm % 
for i=1:dimGATE(1)
    centralDoseGATE = squeeze(dataGATE(i,round(dimGATE(2)/2)-numberVoxels:round(dimGATE(2)/2)+numberVoxels,round(dimGATE(3)/2)-numberVoxels:round(dimGATE(3)/2)+numberVoxels))  ;  
    inDepthDoseGATE(i) = mean2(centralDoseGATE); 
    depthsGATE(i) = (dimGATE(1)-i) .* voxelSize ;
end

depthsGATE = flipud(depthsGATE.');
inDepthDoseGATE = flipud(inDepthDoseGATE.') ;
IDD = flipud(sum(sum(dataGATE,3),2)) ;

figure;
% plot(depthsGATE, inDepthDoseGATE,'-.r' ); hold on; 
plot(depthsGATE, IDD,'-.k' );
% xlim([0 180])
% xticks([0:20:180])
% ylim([0 2])
xlabel('Penentration depth [mm]') 
ylabel('Dose[Gy]')

% get range of proton beam 
Range = getRange(depthsGATE, IDD) ;
Range20 = getRange20(depthsGATE, IDD) ;


% checking lateral profiles at three penetration depths

latCoordinatez = 1:1:dimGATE(3) ;
latCoordinatez = latCoordinatez .* voxelSize - dimGATE(3)./2 - voxelSize./2 ;
d1 = 2 ; % in mm
spotDoseGATE = sum(dataGATE,3) ;
totalDoseSpot = sum(sum(dataGATE,3),2) ;

lateralzGATE_d1 = squeeze(dataGATE(dimGATE(1)-round(d1./voxelSize),round(dimGATE(2)./2),:))  ;
lateralzGATE_d11 = squeeze(dataGATE(dimGATE(1)-round(d1./voxelSize),round(dimGATE(2)./2)+1,:))  ;
lateralzGATE_d12 = squeeze(dataGATE(dimGATE(1)-round(d1./voxelSize),round(dimGATE(2)./2)-1,:))  ;
figure;
plot(latCoordinatez,lateralzGATE_d1) ; hold on; 
plot(latCoordinatez,lateralzGATE_d11) ; hold on; 
plot(latCoordinatez,lateralzGATE_d12) ; hold on; 

            % % fitting lateral profiles
            % % multiple Gaussian 
            % fz1 = fit(latCoordinatez.',lateralzGATE_d1,'gauss1');
            % fz2 = fit(latCoordinatez.',lateralzGATE_d1,'gauss2');
            % fz3 = fit(latCoordinatez.',lateralzGATE_d1,'gauss3');
            % 
            % % 
            % figure;
            % plot(latCoordinatez,lateralzGATE_d1,'+k'); hold on;
            % plot(fz1,'-r'); hold on;
            % plot(fz2,'-b'); hold on;
            % plot(fz3,'-g'); hold on;
            % xlim([-50 50]);
            % ylim([1 1e6]);
            % set(gca,'YScale','log')
            % xlabel('Lateral coordinate [mm]') 
            % ylabel('Dose[Gy]') 

% figure;
% plot(depthsTPS, centralSOBPTPS); hold on;
% plot(depthsGATE, centralSOBPGATE); 

figure;
contourf(centraldoseGATE, 'ShowText', 'off') ;
% xlim([50 150])
% xticks([50:25:150])
% ylim([150 200])
xlabel('Lat. coord. [voxel]') 
ylabel('Depth [voxel]') 
colormap(gca,'jet'); 
% caxis([0 2]);
c2 = colorbar('northoutside');
c2.Label.String = 'Dose (Gy)';
title('GATE'); 

% rVectors = rVectors(:); 
% z = z(:); 
% dose = dose(:);
doseThreshold = 1e-3;

% define the double Gaussian function

fF_Gauss = @(a,x) exp(-(x./(a*2)).^2)/(sqrt(2*pi)*a);
lateralFitFunctionTotal = @(a,x) a(1)*(a(2)*fF_Gauss(a(3),x) + (1-a(2))*fF_Gauss(a(4),x));

% fF_Gauss = @(a,x) exp(-((x./a).^2)/2)./(sqrt(2*pi)*a);
% doubleGaussian = @(a,x) a(1)*(a(2)*fF_Gauss(a(3),x) + (1-a(2))*fF_Gauss(a(4),x));

dataArray = zeros(dimGATE(1), 2);
chiSq = zeros(dimGATE(1), 1);
startPoints = [1,0.9,7,12];


f2GaussParArray = zeros(dimGATE(1),6) ; 
doFigure = false;

% single Gauss fitting parameters for all depths
for i=0:dimGATE(1)-1;  
    depth = i .* voxelSize ; % in mm
     % get parameters only for pen. depths before the R20
            if i < round((Range20)./ voxelSize)
            % check maximum position 
            SpotDoseinDepth = squeeze(dataGATE(dimGATE(1)-i,:,:)) ;
            dataArray(i+1,1) = max(SpotDoseinDepth(:)) ;
            dataArray(i+1,2) = sum(SpotDoseinDepth(:)) ;
            
            normalizedSpotDose = SpotDoseinDepth ./ max(SpotDoseinDepth(:)) ;
            [mxv,idx] = max(normalizedSpotDose(:)); 
            [indexMaxX,indexMaxY,indexMaxZ] = ind2sub(size(normalizedSpotDose),idx) ;
            % lateral profile   
            profX = squeeze(normalizedSpotDose(indexMaxY,:)) ;
            laterallyIntegratedDose = sum(profX);
            normalizedLateralProf = profX ./ laterallyIntegratedDose ;
            
% %                           % restrict fitting to lateral distances where dose is lower than 0.1% of Dmax
                            %%  only consider dose > doseThreshold of total dose
                                [lateralThreashold,lateralDose] = applyDoseThreshold(latCoordinatez,profX,doseThreshold);

                                d = lateralDose ;
                                TotalIntegratedDose = sum(d) ;
%                                 d = d ./ TotalIntegratedDose ;
                                r = lateralThreashold;

                                % get initial parameters from single Gaussian fit
                                  f1Gauss = fit(r.',(d).','gauss1'); 
                                  f1GaussPar = coeffvalues(f1Gauss) ;
                                  f1GaussParArray(i+1,:) = f1GaussPar ;
                
                  %% fit all combinationsOfFunctions - do various Start Points - Andreas software
                    if i > 1
                        startPoints = param(i-1,:,:);
                    else
                        startPoints = [f1GaussPar(1), 0.9, f1GaussPar(3)./2, 1.5 * f1GaussPar(3)./2];
                    end                  
%                      tdoubleGaussian=@(a) lateralFitFunctionTotal(a,r);
%                      [param(i+1,:), chiSq(i+1)] = fitAnyFunction(tdoubleGaussian,d,startPoints);
%                      [~,bestFitInd] = min(chiSq(i+1,:));
%                      lastParValues = param(i+1,:) ;
                     
                tdoubleGaussian=@(a) lateralFitFunctionTotal(a,r);     
                [paramDoseWeighted, ~] = fitAnyFunction(tdoubleGaussian,d,startPoints,true);
                [paramDG, ~] = fitAnyFunction(tdoubleGaussian,d,paramDoseWeighted,true);
                param(i+1,:) = paramDoseWeighted;
                tFixedDoubleGaussian = @(a) a(1)*( a(2)*fF_Gauss(a(3),r) + (1-a(2))*fF_Gauss(paramDoseWeighted(4),r));
                [paramFixedDoubleGaussian, chiSq(i+1)] = fitAnyFunction(tFixedDoubleGaussian,d,paramDoseWeighted(1:4),false);
                param(i+1,:) =[paramFixedDoubleGaussian(1:2), paramDoseWeighted(3:4)];
                currParam = [paramFixedDoubleGaussian, paramDoseWeighted(3)] ;
                [~,bestFitInd] = min(chiSq(i+1,:));
                lastParValues = param(i+1,:) ;
            
            % plotting 
            if doFigure
            figure
            %% plot GATE dose
            semilogy(r, d ,'xk')
            hold on
            %% fit 
            semilogy(r, tdoubleGaussian(param(i+1,:)))
            semilogy(r, tdoubleGaussian(paramDG) )
            semilogy(r, tFixedDoubleGaussian(paramFixedDoubleGaussian) ,'--r' )

%             ylim([doseThreshold,1])

        %     xlim([0,30])
            ylabel('Edep_i / Edep_{max}(x)')
            xlabel('r / mm')
            legend( 'H', 'H fit', 'total', 'total fit')
            title(['Dose at z = ', num2str(i+1), ' mm'])

            hold off
            end
                
            else 
            param(i+1,:) = lastParValues ;  
            end
            

end


for testIndex = 0:dimGATE(1)-1;
% check test at any penetration depth 
  Amp= param(testIndex+1,1) ;
  ratio = param(testIndex+1,2) ;
  Sigma1 = param(testIndex+1,3) ;
  Sigma2 = param(testIndex+1,4) ;
        
  fittedData = Amp*(ratio*fF_Gauss(Sigma1,r) + (1-ratio)*fF_Gauss(Sigma2,r));                  

                    
            % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
            [xi,yi]     = meshgrid(r,r);
            [pdfx,pdfy] = meshgrid(fittedData,fittedData);
            % Calculate combined pdf, under assumption of independence
            pdfxy = pdfx .* pdfy ; 
            norm2Dkernel = pdfxy ./ max(pdfxy(:));
            sum2dpdfKernel = sum(sum(pdfxy,2),1) ;   
                                  % Plot the results
%                                 figure;
%                                 mesh(xi,yi,pdfxy)
%                                 set(gca,'XLim',[min(r) max(r)])
%                                 set(gca,'YLim',[min(r) max(r)])
            
            spotDose = IDD(testIndex+1) * pdfxy   ;                   
            [maxPB,idxPB] = max(spotDose(:)); 
            [indexMaxYPB,indexMaxZPB] = ind2sub(size(spotDose),idxPB) ;
                                
            centraldataGATE = squeeze(dataGATE(dimGATE(1)-testIndex,:,:)) ;
            normdataGATE = centraldataGATE ./ max(centraldataGATE(:)) ;                    
            [maxvalueGATE,idxGATE] = max(centraldataGATE(:)); 
            [indexMaxYGATE,indexMaxZGATE] = ind2sub(size(centraldataGATE),idxGATE) ;                        
                    % check lateral profiles 
                    centralProfYGATE = centraldataGATE(:,indexMaxZGATE);
                    centralProfZGATE = centraldataGATE(indexMaxYGATE,:);
                    centralProfYPB = spotDose(:,indexMaxZPB);
                    centralProfZPB = spotDose(indexMaxYPB,:);
                    
                        depths(testIndex+1) = testIndex .* voxelSize ; % in mm
                        % get parameters only for pen. depths before the
                        % R20
                        if testIndex < round((Range20)./ voxelSize)
                        amplitudeRatioZ(testIndex+1) = max(centralProfZGATE(:)) ./ max(centralProfZPB(:)) ;
                        amplitudeRatioY(testIndex+1) = max(centralProfYGATE(:)) ./ max(centralProfYPB(:)) ;
                        lastamplitudeRatioZ = max(centralProfZGATE(:)) ./ max(centralProfZPB(:)) ;
                        lastamplitudeRatioY = max(centralProfYGATE(:)) ./ max(centralProfYPB(:)) ;
                        else 
                        amplitudeRatioZ(testIndex+1) = lastamplitudeRatioZ ;
                        amplitudeRatioY(testIndex+1) = lastamplitudeRatioY ;
                        end
                    
                           correctedspotDose = amplitudeRatioZ(testIndex+1) .* IDD(testIndex+1) * pdfxy   ;
                           dim2DPB = size(correctedspotDose);
                                    [corrmaxPB,corridxPB] = max(correctedspotDose(:)); 
                                    [corrindexMaxYPB,corrindexMaxZPB] = ind2sub(size(correctedspotDose),corridxPB) ;
                           newcentralProfZPB = correctedspotDose(:,corrindexMaxYPB);
                           newcentralProfYPB = correctedspotDose(corrindexMaxZPB,:);
                           NewamplitudeRatioZ(testIndex+1) = max(centralProfZGATE(:)) ./ max(newcentralProfZPB(:)) ;
                           NewamplitudeRatioY(testIndex+1) = max(centralProfYGATE(:)) ./ max(newcentralProfYPB(:)) ;
                           
                           if testIndex == 2 ;
                              spotPB2DSurf = correctedspotDose ;
                                 [mspotPB2DSurf,idxspotPB2DSurf] = max(spotPB2DSurf(:)); 
                                 [indexMaxYspotPB2DSurf,indexMaxZspotPB2DSurf] = ind2sub(size(spotPB2DSurf),idxspotPB2DSurf) ;
                              spotGATE2DSurf = squeeze(dataGATE(dimGATE(1)-testIndex,:,:)) ;
                                 [mspotGATE2DSurf,idxspotGATE2DSurf] = max(spotGATE2DSurf(:)); 
                                 [indexMaxYspotGATE2DSurf,indexMaxZspotGATE2DSurf] = ind2sub(size(spotGATE2DSurf),idxspotGATE2DSurf) ;
                              centralProfYGATE_test1 = spotGATE2DSurf(:,indexMaxZspotGATE2DSurf);
                              centralProfZGATE_test1 = spotGATE2DSurf(indexMaxYspotGATE2DSurf,:);
                              centralProfYPB_test1 = spotPB2DSurf(:,indexMaxZspotPB2DSurf);
                              centralProfZPB_test1 = spotPB2DSurf(indexMaxYspotPB2DSurf,:);
                              
                           end
                          
                           if testIndex == round(0.5.*Range ./2 ) ;
                              spotPB2D50 = correctedspotDose ;
                                 [mspotPB2D50,idxspotPB2D50] = max(spotPB2D50(:)); 
                                 [indexMaxYspotPB2D50,indexMaxZspotPB2D50] = ind2sub(size(spotPB2D50),idxspotPB2D50) ;
                              spotGATE2D50 = squeeze(dataGATE(dimGATE(1)-testIndex,:,:)) ;
                                 [mspotGATE2D50,idxspotGATE2D50] = max(spotGATE2D50(:)); 
                                 [indexMaxYspotGATE2D50,indexMaxZspotGATE2D50] = ind2sub(size(spotGATE2D50),idxspotGATE2D50) ;
                              centralProfYGATE_test3 = spotGATE2D50(:,indexMaxZspotGATE2D50);
                              centralProfZGATE_test3 = spotGATE2D50(indexMaxYspotGATE2D50,:);
                              centralProfYPB_test3 = spotPB2D50(:,indexMaxZspotPB2D50);
                              centralProfZPB_test3 = spotPB2D50(indexMaxYspotPB2D50,:);
                              
                           end
                           
                           
                          if testIndex == round(0.9.*Range ./2 ) ;
                              spotPB2DBPeak = correctedspotDose ;
                                 [mspotPB2DBPeak,idxspotPB2DBPeak] = max(spotPB2DBPeak(:)); 
                                 [indexMaxYspotPB2DBPeak,indexMaxZspotPB2DBPeak] = ind2sub(size(spotPB2DBPeak),idxspotPB2DBPeak) ;
                              spotGATE2DBPeak = squeeze(dataGATE(dimGATE(1)-testIndex,:,:)) ;
                                 [mspotGATE2DBPeak,idxspotGATE2DBPeak] = max(spotGATE2DBPeak(:)); 
                                 [indexMaxYspotGATE2DBPeak,indexMaxZspotGATE2DBPeak] = ind2sub(size(spotGATE2DBPeak),idxspotGATE2DBPeak) ;
                              centralProfYGATE_test2 = spotGATE2DBPeak(:,indexMaxZspotGATE2DBPeak);
                              centralProfZGATE_test2 = spotGATE2DBPeak(indexMaxYspotGATE2DBPeak,:);
                              centralProfYPB_test2 = spotPB2DBPeak(:,indexMaxZspotPB2DBPeak);
                              centralProfZPB_test2 = spotPB2DBPeak(indexMaxYspotPB2DBPeak,:);
                              
                           end
                           
                           
                           centralDosePB(testIndex+1) = correctedspotDose(round(dim2DPB(1)/2),round(dim2DPB(2)/2))  ;

end ;                                   


figure;
subplot(1,3,1)
contourf(spotPB2DSurf, 'ShowText', 'off','LineStyle','none') ; 
% xlim([80 120])
% xticks([50:25:150])
% ylim([80 120])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,'jet');
% caxis([0 1.1]);
c1 = colorbar('northoutside');
c1.Label.String = 'Dose (Gy)';
title('parametrization - TPS'); 

subplot(1,3,2)
contourf(spotGATE2DSurf, 'ShowText', 'off','LineStyle','none') ;
% xlim([80 120])
% xticks([50:25:150])
% ylim([80 120])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,'jet');
% caxis([0 1.1]);
c2 = colorbar('northoutside');
c2.Label.String = 'Dose (Gy)';
title('GATE'); 


% try fancy colormap ;-)
% Create 13-level blue/red diverging color map:
level = 10; n = ceil(level/2);
cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
cmapBlueRed = [cmap1; cmap2(2:end, :)];
subplot(1,3,3)
spotPB2DSurf(spotPB2DSurf<0.01.*max(spotPB2DSurf(:))) = 0 ;
spotGATE2DSurf(spotGATE2DSurf<0.01.*max(spotGATE2DSurf(:))) = 0 ;
doseDiff = 100 .* (((spotPB2DSurf-spotGATE2DSurf)./spotGATE2DSurf)) ;
% doseDiff(doseDiff == 1) = 0; 
% doseDiff(doseDiff < -1) = -1; 
contourf(doseDiff, 'ShowText', 'off','LineStyle','none') ; 
% xlim([80 120])
% xticks([50:25:150])
% ylim([80 120])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,cmapBlueRed); 
% c1 = min(doseDiff(:));
% c2 = max(doseDiff(:));
caxis([-5 5]);
c3 = colorbar('northoutside');
c3.Label.String = 'dose diff. (%)';
title('abs. dose diff (GATE-TPS)'); 
            
figure;
subplot(1,3,1)
yyaxis left
plot(latCoordinatez,centralProfZPB_test1,'-r') ; hold on;
plot(latCoordinatez,centralProfZGATE_test1,'--b') ;
% xlim([-40 40]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez, 100.*(centralProfZPB_test1-centralProfZGATE_test1)./centralProfZGATE_test1) ; hold on;
ylim([-10 10]);
subplot(1,3,2)
yyaxis left
plot(latCoordinatez,centralProfZPB_test2,'-r') ; hold on;
plot(latCoordinatez,centralProfZGATE_test2,'--b') ;
% xlim([-40 40]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez, 100.*(centralProfZPB_test2-centralProfZGATE_test2)./centralProfZGATE_test2) ; hold on;
ylim([-10 10]);
subplot(1,3,3)
yyaxis left
plot(latCoordinatez,centralProfZPB_test3,'-r') ; hold on;
plot(latCoordinatez,centralProfZGATE_test3,'--b') ;
% xlim([-40 40]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez, 100.*(centralProfZPB_test3-centralProfZGATE_test3)./centralProfZGATE_test3) ; hold on;
ylim([-10 10]);


figure;
subplot(1,3,1)
contourf(spotPB2DBPeak, 'ShowText', 'off','LineStyle','none') ; 
% xlim([80 120])
% xticks([50:25:150])
% ylim([80 120])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,'jet');
% caxis([0 1.1]);
c1 = colorbar('northoutside');
c1.Label.String = 'Dose (Gy)';
title('parametrization - TPS'); 

subplot(1,3,2)
contourf(spotGATE2DBPeak, 'ShowText', 'off','LineStyle','none') ;
% xlim([80 120])
% xticks([50:25:150])
% ylim([80 120])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,'jet');
% caxis([-0.1 0.1]);
c2 = colorbar('northoutside');
c2.Label.String = 'Dose (Gy)';
title('GATE'); 


% try fancy colormap ;-)
% Create 13-level blue/red diverging color map:
level = 10; n = ceil(level/2);
cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
cmapBlueRed = [cmap1; cmap2(2:end, :)];

subplot(1,3,3)
potPB2DBPeak(spotPB2DBPeak<0.01.*max(spotPB2DBPeak(:))) = 0 ;
spotGATE2DBPeak(spotGATE2DBPeak<0.01.*max(spotGATE2DBPeak(:))) = 0 ;
doseDiff2 = 100 .* (((spotPB2DBPeak-spotGATE2DBPeak)./spotGATE2DBPeak)) ;
% doseDiff(doseDiff == 1) = 0; 
% doseDiff(doseDiff < -1) = -1; 
contourf(doseDiff2, 'ShowText', 'off','LineStyle','none') ; 
% xlim([70 120])
% xticks([50:25:150])
% ylim([70 120])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,cmapBlueRed); 
caxis([-100 100]);
c3 = colorbar('northoutside');
c3.Label.String = 'dose diff. (Gy)';
title('abs. dose diff (GATE-TPS)'); 
            


%% create ascii output file %%
for i=1:size(param,1)
AmpArray = param(i,1) ;
ratioArray = param(i,2) ;
Sigma1Array = param(i,3) ;
Sigma2Array = param(i,3) ; 

end

outputDataZ = [depthsGATE IDD param amplitudeRatioZ.'] ;
dlmwrite('ParametersDoubleGauss_IR1BioSetup-100.4MeV.txt',outputDataZ,'delimiter','\t') ;


%% check gamma index pass-rates %%
ctRes = [1 1] ;
gammaIndexCriteria = [2 2]; 
% testIndex1 = 20 ;
interpPoints = 0 ;
slice = 1 ;

outputData1 = matRad_gammaIndex2D(spotPB2DSurf,spotGATE2DSurf,ctRes,slice,gammaIndexCriteria,interpPoints,'global') ;
outputData2 = matRad_gammaIndex2D(spotPB2D50,spotGATE2D50,ctRes,slice,gammaIndexCriteria,interpPoints,'global') ;
outputData3 = matRad_gammaIndex2D(spotPB2DBPeak,spotGATE2DBPeak,ctRes,slice,gammaIndexCriteria,interpPoints,'global') ;


function [r_forDoseAboveThreshold,doseAboveThreshold] = applyDoseThreshold(r,d,dThreshold)
    dM = d>dThreshold;
    r_forDoseAboveThreshold = r(dM);
    doseAboveThreshold = d(dM);
end

function doseInRNormalized = getDoseInZnormalized(d,z,dTotal)
    doseInR = squeeze(d(:,z));
    doseInRTotal = squeeze(dTotal(:,z));
    doseInRNormalized = doseInR./max(doseInRTotal);
end


