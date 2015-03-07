function ienergy = getienergytest(pp,gr,yi,na)

yiext = [yi 0 ];
yiprime = diff(yiext)./gr.dx;

ienergy = na*2*pi*sum(gr.dx.*gr.x.*((1/2)*pp.alpha*abs(yiprime).^2+(1/2).*(1/pp.alpha).*(pp.omegat^2).*gr.x.^2.*(yi.^2)+(1/2)*pp.etaa*(na-1)*pp.beta*abs(yi).^4));

end

