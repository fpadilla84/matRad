function [param, chiSq] = fitAnyFunction(functHandl,y, startPoints,doOneOverAbsValueFlag)
if nargin <4
    doOneOverAbsValueFlag = false;
end
    if doOneOverAbsValueFlag
        myFitFunct = @(a) ((functHandl(a) - y)./((y)));
    else
        myFitFunct = @(a) ((functHandl(a) - y));
    end
    numelParam = numel(startPoints);
    lb = zeros(1,numelParam);
    ub = 100*ones(1,numelParam);
    options = optimoptions('lsqnonlin','MaxFunEvals', 50,'MaxIter',50, 'TolFun', 1e-11,'StepTolerance',1e-11,'Display','off');
    if numelParam >2
        ub(2) = 0.9999;
    end
    if numelParam > 4
        ub(2:3) = 0.999;
    end
    [param, chiSq] = lsqnonlin(myFitFunct, startPoints ,lb,ub,options);
%     chiSqT = sum(myFitFunct(param).^2);
end
