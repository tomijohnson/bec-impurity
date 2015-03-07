function [sol,E_min] = minimizetp(param,sol0)

    a = param(1);
    b = param(2);
    d = param(3);
    r = param(4);
    N = param(5);

    fun = @(sol) E(sol);
    
    [sol,E_min] = fminsearch(fun,sol0); 

    function Energy = E(sol)
        Energy = a*(1+(N-1)*r)/sol(1)^2+(d/N)*(3*sol(2)^2-log(sol(2)^2))+b/(1+(sol(2)/sol(1))^2);
    end

end