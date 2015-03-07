function [E_min,sol] = energyvsN(param,sol0)
        
    E_min = zeros(1,3);
    sol = zeros(3,2);
    for N = 1:3
        [solN,E_minN] = minimizetp([param N],sol0);
        E_min(N) = E_minN;
        sol(N,:) = solN;
    end
    
end