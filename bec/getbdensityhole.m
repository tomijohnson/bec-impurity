function vdensityhole = getbdensityhole(pp,gr,yb)
vdensityhole = 2*pi*sum(gr.dx.*gr.x.*(1/pp.beta - yb.^2));
end