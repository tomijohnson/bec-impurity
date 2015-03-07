function vgtpenergy = getvgenergytp1(pp,na,a,w)
    vgtpenergy = na*pp.etaab*(1+((w/a)^2)*exp((w/a)^2)*(-expint((w/a)^2)));
end