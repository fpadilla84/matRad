function doubleGauss = doubleGaussFunction(x,A,Mean1,Sigma1,w,Mean2,Sigma2)
% %% fpadilla: double Gauss function 
% TAILEDGAUSS produces gauss with exponential tail
%   used in the non-linear fitting
% Parameters 
% Mean1,2          - position of the Gaussian peak2
% Sigma1,2         - width of the Gaussian2
% A                - Normalization Areas corresponding to each Gaussian
% ratio            - Area fraction weights 
pdfSG1    =  @(x,Mean1,Sigma1) (1/2.*pi/(Sigma1^2)).*exp(-(x-Mean1).^2./(Sigma1^2));
pdfSG2    =  @(x,Mean2,Sigma2) (1/2.*pi/(Sigma2^2)).*exp(-(x-Mean2).^2./(Sigma2^2));

% Gauss1    =  @(x,Mean1,Sigma1) exp(- 0.5 .* (x - Mean1).^2 ./ (Sigma1.^2)) ./ sqrt(2.*pi*Sigma1.^2) ;
% Gauss2    =  @(x,Mean1,Sigma1) exp(- 0.5 .* (x - Mean2).^2 ./ (Sigma2.^2)) ./ sqrt(2.*pi*Sigma2.^2) ;
% Gauss2    =  @(x,Mean2,Sigma2) exp(-(x - Mean2).^2 ./ (Sigma2.^2)) ./ sqrt(pi*Sigma1.^2) ;
% ratio = Amp2 ./ Amp1   ;
% Combining the three components

doubleGauss = A .* ((1-w)* pdfSG1(x,Mean1,Sigma1) + w * pdfSG2(x,Mean2,Sigma2));

end

