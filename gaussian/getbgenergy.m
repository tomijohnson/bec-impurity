function vgenergy = getbgenergy(pp,gr,yb,na,a)

vgenergy = na*pp.etaab*pp.beta*2*pi*sum(gr.dx.*gr.x.*getgdensity(gr.x,a).*(abs(yb).^2));

end