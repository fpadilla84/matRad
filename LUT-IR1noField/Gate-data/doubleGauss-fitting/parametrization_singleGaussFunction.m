function singleGauss = parametrization_singleGaussFunction(x,H,Mean,Sigma)
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

pdfSG    =  @(x,Mean,Sigma) (1/sqrt(pi)/Sigma).*exp(-(x-Mean).^2./(Sigma^2));

% Gauss2    =  @(x,Mean2,Sigma2) exp(-(x - Mean2).^2 ./ (Sigma2.^2)) ./ sqrt(pi*Sigma2.^2) ;
% Gauss2    =  @(x,Mean2,Sigma2) exp(-(x - Mean2).^2 ./ (Sigma2.^2)) ./ sqrt(pi*Sigma1.^2) ;
% ratio = Amp2 ./Amp1   ;
% Combining the three components

singleGauss = H * pdfSG(x,Mean,Sigma) ;

end

