function gdensity=getgdensity(x,a)
    %BC: Evaluates the residue of the boundary condition
    gdensity=(1/(pi*a^2))*exp(-(x/a).^2);
end