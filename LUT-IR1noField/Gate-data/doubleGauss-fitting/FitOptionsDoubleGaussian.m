%% Create Default Fit Options for Gaussian Fit  

% Copyright 2015 The MathWorks, Inc.


%%  
options = fitoptions('gauss2')   

  Normalize: 'off'
          Exclude: []
          Weights: []
           Method: 'NonlinearLeastSquares'
           Robust: 'Off'
       StartPoint: [1x0 double]
            Lower: [-Inf -Inf 0 -Inf -Inf 0]
            Upper: [1x0 double]
        Algorithm: 'Trust-Region'
    DiffMinChange: 1.0000e-08
    DiffMaxChange: 0.1000
          Display: 'Notify'
      MaxFunEvals: 600
          MaxIter: 400
           TolFun: 1.0000e-06
             TolX: 1.0000e-06
