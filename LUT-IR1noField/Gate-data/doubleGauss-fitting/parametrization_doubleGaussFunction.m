function doubleGauss = parametrization_doubleGaussFunction(x,H1,Mean1,Sigma1,H2,Mean2,Sigma2)
% %% fpadilla: new expGauss function 
% TAILEDGAUSS produces gauss with exponential tail
%   used in the non-linear fitting
% Parameters 
% Mean          - position of the primary Gaussian peak
% Sigma         - width of the Gaussian
% LGamma        - width of the left exponential
% RGamma        - width of the right exponential
% Lw            - Area fraction corresponding to the left exponential tail
% Rw            - Area fraction corresponding to the right exponential tail

pdfSG1    =  @(x,Mean1,Sigma1) (1/sqrt(pi)/Sigma1).*exp(-(x-Mean1).^2./(Sigma1^2));
pdfSG2    =  @(x,Mean2,Sigma2) (1/sqrt(pi)/Sigma2).*exp(-(x-Mean2).^2./(Sigma2^2));

% Gauss1    =  @(x,Mean1,Sigma1) exp(- 0.5 .* (x - Mean1).^2 ./ (Sigma1.^2)) ./ sqrt(2.*pi*Sigma1.^2) ;
% Gauss2    =  @(x,Mean1,Sigma1) exp(- 0.5 .* (x - Mean2).^2 ./ (Sigma2.^2)) ./ sqrt(2.*pi*Sigma2.^2) ;
% Gauss2    =  @(x,Mean2,Sigma2) exp(-(x - Mean2).^2 ./ (Sigma2.^2)) ./ sqrt(pi*Sigma1.^2) ;
% ratio = Amp2 ./ Amp1   ;
% Combining the three components

doubleGauss = H1 * pdfSG1(x,Mean1,Sigma1) + H2 * pdfSG2(x,Mean2,Sigma2);

end

