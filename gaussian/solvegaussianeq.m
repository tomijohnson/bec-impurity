function  a=solvegaussianeq(pp,np,gr,yb,na,ainit)
    fun = @f; % impurity euler lagrange equation is f = 0
    itcount = 0;
    iterate = 1;
    ainit = ainit*[1 linspace(1/5,5,np.maxit2-1)];
    while iterate == 1
        itcount = itcount + 1;
        a = fzero(fun,ainit(itcount));
        
        if (100 > a) && (a > 1e-14)
            iterate = 0;
        end
        
        if itcount >= np.maxit2
            iterate = 0;
        end
    end
    function yg = f(a)
        yg = (1/pp.alpha)*(pp.omegat^2)*a^4 -(pp.alpha + pp.etaa*(na-1)*pp.beta*(1/(2*pi))) + 4*pp.etaab*pp.beta*sum(gr.dx.*(gr.x).*((gr.x/a).^2-1).*exp(-(gr.x/a).^2).*(abs(yb).^2));
    end
end
