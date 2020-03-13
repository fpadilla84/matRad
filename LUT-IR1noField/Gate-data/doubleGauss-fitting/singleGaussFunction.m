function pdfGauss = singleGaussFunction(x,A,Mean,Sigma)
% %% fpadilla: single Gauss function 
%   used in the non-linear fitting
% Parameters 
% Mean          - position of the primary Gaussian peak
% Sigma         - width of the Gaussian
% Rw            - Total normalization area 

pdfSG1    =  @(x,A,Mean,Sigma) (1/(2.*pi)./(Sigma.^2)).*exp(-(x-Mean).^2./(Sigma^2));

pdfGauss = pdfSG1(x,A,Mean,Sigma) ;

end

