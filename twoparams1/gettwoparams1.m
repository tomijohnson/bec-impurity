function  rmin=gettwoparams1(pp,na,rinit)
    fun = @f;
    
    rmin = fminsearch(fun,rinit);
    
    function fval = f(r)
        fval = na*((pp.alpha/2)+(na-1)*pp.etaa*(pp.beta/2)*(1/(2*pi)))/r(1)^2+(pi/(2*pp.beta))*(r(2)^2-log(r(2)^2))+pp.etaab*na*(1+((r(2)/r(1))^2)*exp((r(2)/r(1))^2)*(-expint((r(2)/r(1))^2)));
    end
end