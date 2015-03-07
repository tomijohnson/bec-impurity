function vienergy = getbienergy(pp,gr,yb,yi,na)

vienergy = na*pp.etaab*pp.beta*2*pi*sum(gr.dx.*gr.x.*(abs(yi).^2).*(abs(yb).^2));
end