close all
clear all

[~,userName]=system('whoami');
if ~isempty(regexp(userName,'fpadilla','Once'))
    tPathGATE='/home/fpadilla/Simulations/IR1-BeamModel/BioSetup/72.4MeV';

else ~isempty(regexp(userName,'fpc','Once'))
     tPathGATE='/Users/fpc/Documents/GATE-Simulations/BioSetup/72.4MeV';
    
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

% define the double Gaussian function
singleGauss = @(a,x) exp(-((x./a).^2)/2)./(sqrt(2*pi)*a);
doubleGaussian = @(a,x) a(1)*( a(2)*singleGauss(a(3),x) + (1-a(2))*singleGauss(a(4),x));

f2GaussParArray = zeros(dimGATE(1),6) ; 

figure;
% single Gauss fitting parameters for all depths
for i=0:dimGATE(1)-1;  
    d = i .* voxelSize ; % in mm
    % get parameters only for pen. depths before the beam range + 5 mm
            if i < round(Range20./ voxelSize)
            % z coordinate - where deflection is observed    
            lateralzGATE = squeeze(dataGATE(dimGATE(1)-i,round(dimGATE(2)./2),:)) .' ;
            % restrict fitting to lateral distances where dose is lower
            % than 0.1% of Dmax
            outliers = excludedata(latCoordinatez,lateralzGATE,'range',[0.95*max(lateralzGATE(:)) max(lateralzGATE(:))]);
            % get initial parameters from single Gaussian fit
            f1Gauss = fit(latCoordinatez.',lateralzGATE.','gauss1','Exclude',outliers,'Robust','on'); 
            f1GaussPar = coeffvalues(f1Gauss) ;
                              plot(latCoordinatez,lateralzGATE,'+k'); hold on;
                              plot(f1Gauss,'-r'); hold on;
             % inital guess parameters 
             
%             f2GaussParArray(i+1,:) = f1GaussPar ;
%             lastParValues = f1GaussPar ;
            % add a second Gaussian 
%             f2Gauss = fit(latCoordinatez.',lateralzGATE.','gauss2');
                      
            
            optionsLatz = fitoptions('gauss2', 'StartPoint',[f1GaussPar(1) f1GaussPar(2) f1GaussPar(3) 0.1.*f1GaussPar(1) f1GaussPar(2) 2.*f1GaussPar(3)], ...   
                                      'Lower', [0.8.*f1GaussPar(1) f1GaussPar(2) f1GaussPar(3) 0.2.*f1GaussPar(1) f1GaussPar(2) f1GaussPar(3)], ...
                                      'Upper', [f1GaussPar(1) f1GaussPar(2) f1GaussPar(3) 0.5.*f1GaussPar(1) f1GaussPar(2) 2.*f1GaussPar(3)]) ; 
            f2Gauss = fit(latCoordinatez.',lateralzGATE.','gauss2',optionsLatz);
%             f2Gauss = fit(latCoordinatez.',lateralzGATE.','gauss2');
%             %%% check parameters 
            f2GaussPar = coeffvalues(f2Gauss) ; 
            f2GaussParArray(i+1,:) = f2GaussPar ;
            lastParValues = f2GaussPar ;
%           

                  plot(f2Gauss,'-b'); hold on; 

            else 
            f2GaussParArray(i+1,:) = lastParValues ;  
            end
    
end

% check test at any penetration depth 
testIndex = 10 ; 
lateralzGATE_testIndex = squeeze(dataGATE(dimGATE(1)-testIndex,round(dimGATE(2)./2),:)) .' ;
intzGate_testIndex = sum(lateralzGATE_testIndex) ;
normlateralzGATE_testIndex = lateralzGATE_testIndex ./ max(lateralzGATE_testIndex(:))  ;
zAmp1= f2GaussParArray(testIndex,1) ;
zMean1 = f2GaussParArray(testIndex,2) ;
zSigma1 = f2GaussParArray(testIndex,3) ;
zH1 = zAmp1 * zSigma1 * sqrt(pi) ;
zAmp2 = f2GaussParArray(testIndex,4) ;
zMean2 = f2GaussParArray(testIndex,5) ;
zSigma2 = f2GaussParArray(testIndex,6) ;
zH2 = zAmp2 * zSigma2 * sqrt(pi) ;
fittedDataz = parametrization_doubleGaussFunction(latCoordinatez, zH1,zMean1,zSigma1,zH2,zMean2,zSigma2) ;
intFittedDataz_testIndex = sum(fittedDataz) ;
normfittedDataz = fittedDataz ./ max(fittedDataz(:)) ;
ratioz = intFittedDataz_testIndex ./ intzGate_testIndex ;

lateralyGATE_testIndex = squeeze(dataGATE(dimGATE(1)-testIndex,:,round(dimGATE(3)./2)))  ;
normlateralyGATE_testIndex = lateralyGATE_testIndex ./ max(lateralyGATE_testIndex(:)) ;
yAmp1= f2GaussParArray(testIndex,1) ;
yMean1 = f2GaussParArray(testIndex,2) ;
ySigma1 = f2GaussParArray(testIndex,3) ;
yH1 = yAmp1 * ySigma1 * sqrt(pi) ;
yAmp2 = f2GaussParArray(testIndex,4) ;
yMean2 = f2GaussParArray(testIndex,5) ;
ySigma2 = f2GaussParArray(testIndex,6) ;
yH2 = yAmp2 * ySigma2 * sqrt(pi) ;
fittedDatay = parametrization_doubleGaussFunction(latCoordinatez, yH1,yMean1,ySigma1,yH2,yMean2,ySigma2) ;
normfittedDatay = fittedDatay ./ max(fittedDatay(:)) ;

ratio = yAmp2 ./ yAmp1 ;

figure;
subplot(1,2,1)
yyaxis left
plot(latCoordinatez,fittedDataz,'-r') ; hold on;
plot(latCoordinatez,lateralzGATE_testIndex,'--b') ;
xlim([-30 30]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez,100.*(fittedDataz-lateralzGATE_testIndex)./lateralzGATE_testIndex) ; hold on;
ylim([-10 10]);

subplot(1,2,2)
yyaxis left
plot(latCoordinatez,fittedDatay,'-r') ; hold on;
plot(latCoordinatez,lateralyGATE_testIndex,'--b') ;
xlim([-30 30]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez,100.*(fittedDatay-lateralyGATE_testIndex)./ lateralyGATE_testIndex) ; hold on;
ylim([-10 10]);


% % check 2D distribution 
% xi = -199:2:199 ; 
% yi = -199:2:199 ;
% pdfxG1 = parametrization_singleGaussFunction(xi, zMean1,zSigma1);
% pdfyG1 = parametrization_singleGaussFunction(yi, yMean1,ySigma1);
% pdfxG2 = parametrization_singleGaussFunction(xi, zMean2,zSigma2);
% pdfyG2 = parametrization_singleGaussFunction(yi, yMean2,ySigma2);
% % Create 2-d grid of coordinates and function values, suitable for 3-d plotting
%             [xxi,yyi]     = meshgrid(xi,yi);
%             [pdfxxG1,pdfyyG1] = meshgrid(pdfxG1,pdfyG1);
%             [pdfxxG2,pdfyyG2] = meshgrid(pdfxG2,pdfyG2);
%             % Calculate combined pdf, under assumption of independence
%             pdfxy = (1-ratio) .* pdfxxG1 .* pdfyyG1 + ratio .* pdfxxG2 .* pdfyyG2  ; 
%             spotDose = IDD(testIndex) .* pdfxy ;
%             norm2Dkernel = pdfxy ./ max(pdfxy(:));
%             sum2dpdfKernel = sum(sum(pdfxy,2),1) ;
%             % Plot the results
%             figure;
%             mesh(xxi,yyi,pdfxy)
%             set(gca,'XLim',[min(xi) max(xi)])
%             set(gca,'YLim',[min(yi) max(yi)])
            
% check 2D distribution 
zi = -199:2:199 ; 
zcentral= -3*zSigma1:2:3.*zSigma1 ;
yi = -199:2:199 ;
ycentral= -3*ySigma1:2:3.*ySigma1 ;


% pdfy = parametrization_doubleGaussFunction(latCoordinatez, yH1,yMean1,ySigma1,yH2,yMean2,ySigma2) ;
% pdfz = parametrization_doubleGaussFunction(latCoordinatez, zH1,zMean1,zSigma1,zH2,zMean2,zSigma2) ;
% totalAreapdfz = sum(pdfy) ;
% totalAreapdfy = sum(pdfz) ;
% 
% % 
% pdfy = pdfy ./ totalAreapdfz  ;
% pdfz = pdfz ./ totalAreapdfy ;

zH = zH1 ;
ratioz = zSigma2 .* zAmp2 ./ zSigma1 ./ zAmp1 ;
yH = yH1 ;
ratioy = ySigma2 .* yAmp2 ./ ySigma1 ./ yAmp1 ;
% 
pdfz = doubleGaussFunction(zi, zMean1,zSigma1,ratioz,zMean2,zSigma2);
pdfy = doubleGaussFunction(yi, yMean1,ySigma1,ratioy,yMean2,ySigma2);
% 
pdfz1 = singleGaussFunction(zi, zMean1,zSigma1) ;
pdfy1 = singleGaussFunction(yi, yMean1,ySigma1) ; 
pdfz2 = singleGaussFunction(zi, zMean2,zSigma2) ;
pdfy2 = singleGaussFunction(yi, yMean2,ySigma2) ; 

alphaz = zSigma1 .* zAmp1 ./ zSigma2 ./ zAmp2 ;
ratioz = 1 ./ (alphaz + 1) ;
alphay = ySigma1 .* yAmp1 ./ ySigma2 ./ yAmp2 ;
ratioy = 1 ./ (alphay + 1) ;

pdfz = (1-ratioz) .* pdfz1 + ratioz .* pdfz2 ;
pdfy = (1-ratioy) .* pdfy1 + ratioy .* pdfy2 ;
% 
% areaRatioz = sum(pdfz(size(pdfz,2)./2-round(2.*zSigma1./voxelSize):size(pdfz,2)./2+round(2.*zSigma1./voxelSize))) ./ sum(pdfz) ;
% areaRatioy = sum(pdfy(size(pdfy,2)./2-round(2.*ySigma1./voxelSize):size(pdfy,2)./2+round(2.*ySigma1./voxelSize))) ./ sum(pdfy) ;

% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
            [zzi,yyi]     = meshgrid(zi,yi);
            [pdfzz,pdfyy] = meshgrid(pdfz,pdfy);
            % Calculate combined pdf, under assumption of independence
            pdfzy = pdfzz .* pdfyy ; 
            spotDose = 4 .* IDD(testIndex) * pdfzy   ;
            norm2Dkernel = pdfzy ./ max(pdfzy(:));
            sum2dpdfKernel = sum(sum(pdfzy,2),1) ;
            % Plot the results
            figure;
            mesh(zzi,yyi,pdfzy)
            set(gca,'XLim',[min(zi) max(zi)])
            set(gca,'YLim',[min(yi) max(yi)])


centraldataGATE = squeeze(dataGATE(dimGATE(1)-testIndex,:,:)) ;
normdataGATE = centraldataGATE ./ max(centraldataGATE(:)) ;

centraldataGATE(centraldataGATE<0.001.*max(centraldataGATE(:))) = 0;
spotDose(spotDose<0.001.*max(spotDose(:))) = 0;

centralProfZGATE = centraldataGATE(:,size(centraldataGATE,2)./2);
centralProfYGATE = centraldataGATE(size(centraldataGATE,1)./2,:);
centralProfZPB = spotDose(:,size(spotDose,2)./2);
centralProfYPB = spotDose(size(spotDose,1)./2,:);


figure;
subplot(1,3,1)
contourf(norm2Dkernel, 'ShowText', 'off','LineStyle','none') ; 
xlim([90 110])
% xticks([50:25:150])
ylim([90 110])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,'jet');
caxis([0 1.1]);
c1 = colorbar('northoutside');
c1.Label.String = 'Dose (Gy)';
title('parametrization - TPS'); 

subplot(1,3,2)
contourf(normdataGATE, 'ShowText', 'off','LineStyle','none') ;
xlim([90 110])
% xticks([50:25:150])
ylim([90 110])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,'jet');
caxis([0 1.1]);
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
doseDiff = ((norm2Dkernel-normdataGATE)) ;
doseDiff(doseDiff == 1) = 0; 
doseDiff(doseDiff < -1) = -1; 
contourf(doseDiff, 'ShowText', 'off','LineStyle','none') ; 
xlim([90 110])
% xticks([50:25:150])
ylim([90 110])
xlabel('Lat. coord. [voxel]') 
ylabel('Lat. coord. [voxel]') 
colormap(gca,cmapBlueRed); 
caxis([-0.05 0.05]);
c3 = colorbar('northoutside');
c3.Label.String = 'dose diff. (Gy)';
title('abs. dose diff (GATE-TPS)'); 
            
figure;
subplot(1,2,1)
yyaxis left
plot(latCoordinatez,centralProfZPB ./ max(centralProfZPB(:)),'-r') ; hold on;
plot(latCoordinatez,centralProfZGATE ./ max(centralProfZGATE(:)),'--b') ;
xlim([-30 30]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez, 100.*(centralProfZPB-centralProfZGATE)./centralProfZGATE) ; hold on;
ylim([-10 10]);
subplot(1,2,2)
yyaxis left
plot(latCoordinatez,centralProfYPB,'-r') ; hold on;
plot(latCoordinatez,centralProfYGATE,'--b') ;
xlim([-30 30]);
% ylim([0 1.1]);
yyaxis right
plot(latCoordinatez, 100.*(centralProfYPB-centralProfYGATE)./centralProfYGATE) ; hold on;
ylim([-10 10]);


%% create ascii output file %%
for i=1:size(f2GaussParArray,1)
zAmp1Array = f2GaussParArray(i,1) ;
zMean1Array = f2GaussParArray(i,2) ;
zSigma1Array = f2GaussParArray(i,3) ;
zAmp2Array = f2GaussParArray(i,4) ;
zMean2Array = f2GaussParArray(i,5) ;
zSigma2Array = f2GaussParArray(i,6) ;
yAmp1Array = f2GaussParArray(i,1) ;
yMean1Array = f2GaussParArray(i,2) ;
ySigma1Array = f2GaussParArray(i,3) ;
yAmp2Array = f2GaussParArray(i,4) ;
yMean2Array = f2GaussParArray(i,5) ;
ySigma2Array = f2GaussParArray(i,6) ;

alphazArray = zSigma1Array .* zAmp1Array ./ zSigma2Array ./ zAmp2Array ;
alphayArray = ySigma1Array .* yAmp1Array ./ ySigma2Array ./ yAmp2Array ;

ratiozArray(i) = 1 ./ (alphazArray + 1) ;
ratioyArray(i) = 1 ./ (alphayArray + 1) ;
end
outputDataZ = [depthsGATE IDD f2GaussParArray ratiozArray.'] ;
outputDataY = [depthsGATE IDD f2GaussParArray ratioyArray.'] ;

dlmwrite('ParametersZDoubleGauss_216MeV.txt',outputDataZ,'delimiter','\t') ;
dlmwrite('ParametersYDoubleGauss_216MeV.txt',outputDataY,'delimiter','\t') ;

%% check gamma index pass-rates %%
ctRes = [2 2] ;
gammaIndexCriteria = [2 2]; 
slice = 100 ;
interpPoints = 0 ;

outputData = matRad_gammaIndex2D(norm2Dkernel,normdataGATE,ctRes,slice,gammaIndexCriteria,interpPoints,'global') ;


