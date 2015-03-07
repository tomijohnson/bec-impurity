function genergy = getgenergy(pp,na,a)

genergy = na*((1/2)*pp.alpha*1/a^2 + (1/2)*(1/pp.alpha)*(pp.omegat^2)*(a^2) + (1/2)*pp.etaa*(na-1)*pp.beta/(2*pi*a^2));
end