function [estimates0, model0] = fitcurve0(xdata, ydata,start_point)
model0 = @expfun0;
estimates0 = fminsearch(model0, start_point);%,optimset('MaxFunEvals',1000000));
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = expfun0(params)
        A = params(1);
        B = params(2);
        FittedCurve = A*exp(-B*xdata)./(xdata(end).^(1/2));
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector.^ 2);
    end
end