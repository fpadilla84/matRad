function expGauss = matRad_expGaussFunction(vR,Mean,Sigma,LGamma,RGamma,Lw,Rw)
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

Gauss    =  @(vR,Mean,Sigma) exp(-(vR - Mean).^2 ./ (2*Sigma.^2))./ sqrt(2*pi*Sigma.^2) ;
LeftTail = @(vR,Mean,Sigma,LGamma) (0.5./LGamma./Sigma./exp(-0.5./LGamma.^2)).* exp((vR-Mean)./LGamma./Sigma).* ...
    erfc(((vR-Mean)./Sigma./sqrt(2)) + 1./LGamma./sqrt(2)) ;
RightTail = @(vR,Mean,Sigma,RGamma) (0.5./RGamma./Sigma./exp(-0.5./RGamma.^2)).* exp((Mean-vR)./RGamma./Sigma).* ...
    erfc(((Mean-vR)./Sigma./sqrt(2)) + 1./RGamma./sqrt(2)) ;

% Combining the three components

expGauss = (1-Lw-Rw).* Gauss(vR,Mean,Sigma) + Lw .* LeftTail(vR,Mean,Sigma,LGamma) + Rw .* RightTail(vR,Mean,Sigma,RGamma);

end

