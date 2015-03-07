function ienergy = getienergy(pp,gr,yi,na)

yiext = [yi 0 ];
yiprime = diff(yiext)./gr.dx;

yiext = [0 yiprime];
yi2prime = diff(yiext)./gr.dx;

ienergy = na*2*pi*sum(gr.dx.*(-(1/2)*pp.alpha*gr.x.*conj(yi).*yi2prime-(1/2)*pp.alpha*conj(yi).*yiprime+gr.x.*(1/2).*(1/pp.alpha).*(pp.omegat^2).*gr.x.^2.*(yi.^2)+(1/2)*pp.etaa*(na-1)*pp.beta*gr.x.*abs(yi).^4));


end