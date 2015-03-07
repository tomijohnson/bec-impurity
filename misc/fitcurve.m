function [estimates, model] = fitcurve(xdata, ydata,start_point)
model = @expfun;
estimates = fminsearch(model, start_point);%,optimset('MaxFunEvals',1000000));
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = expfun(params)
        A = params(1);
        B = params(2);
        FittedCurve = A*exp(-B*xdata)./(xdata.^(1/2));
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector.^ 2);
    end
end