function [yb,solb]=solvebeceqg(pp,np,gr,na,a)
    x = gr.x;
    l = min(np.maxx,pp.becrad);
    id = find(x<l);
    xred = [x(id) l];
    
    if (pp.vortex == 1)
        %this solves the vortex equation with the gaussian impurity
        solb=bvpinit(xred,@guessvortex,1);
        options = bvpset('RelTol',1e-4,'SingularTerm',[0 0;0 -3]);
        solb=bvp4c(@bvpfunvortex,@bcvortex,solb,options);
        sxint = deval(solb,xred);
        sxint(1,:) = xred.*sxint(1,:);
    else if (pp.vortex == 0)
        %this solves the bec equation with the gaussian impurity
        solb=bvpinit(xred,@guess,1);
        options = bvpset('RelTol',1e-4,'SingularTerm',[0 0;0 -1]);
        solb=bvp4c(@bvpfun,@bc,solb,options);
        sxint = deval(solb,xred);
        sxint(1,:) = sxint(1,:);
        end
    end

%     solv.parameters
    yb = [sxint(1,1:end-1) (1/sqrt(pp.beta))*ones(1,size(x,2)-size(xred,2)+1)];


    function ybprime = bvpfun(x,yb,mub)
        ybprime=[yb(2); 2*(pp.etaab*na*pp.beta*getgdensity(x,a)- mub)*yb(1) + 2*pp.beta*yb(1)^3];
    end

    function ybprime = bvpfunvortex(x,yb,mub)
        ybprime=[yb(2); 2*(pp.etaab*na*pp.beta*getgdensity(x,a)- mub)*yb(1) + 2*pp.beta*(x^2)*yb(1)^3];
    end
    
    function res=bc(yb0,yb1,~)
        res=[yb1(1)-1/sqrt(pp.beta); yb1(2); yb0(2)];
    end

    function res=bcvortex(yb0,yb1,~)
        res=[yb1(1)-1/(sqrt(pp.beta)*l); yb1(2)+1/(sqrt(pp.beta)*l^2); yb0(2)];
    end
    
    function yinit = guess(x) 
        yinit = (1/sqrt(pp.beta))*0.5*[1+tanh(x).^2 2*tanh(x).*sech(x).^2];
    end

    function ybinit = guessvortex(x) 
        ybinit = (1/sqrt(pp.beta))*[1/(1+x^2)^(1/2) (-x/(1+x^2)^(3/2))];
    end
end