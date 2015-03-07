function [estimates, model] = fitcurvevortextp(xdata, ydata,start_point)
model = @expfun;
estimates = fminsearch(model, start_point);%,optimset('MaxFunEvals',1000000));
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = expfun(params)
        A = params(1);
        B = params(2);
        FittedCurve = (1./(1+(A./xdata).^2)).*(1-exp(-(xdata/B).^2));
        ErrorVector = FittedCurve - ydata.^2;
        sse = sum(ErrorVector.^ 2);
    end
end