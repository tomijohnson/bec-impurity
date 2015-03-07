function benergy = getbenergy(pp,gr,yb)

%extract wavefunction and its derivatives
ybext = [yb yb(end)];
ybprime = diff(ybext)./gr.dx;

ybext = [ybprime ybprime(end)];
yb2prime = diff(ybext)./gr.dx;

%calculate energy (1/2)*abs(yvprime).^2 +
benergy = 2*pi*sum(gr.dx.*gr.x.*((yb.^2).*((1/2)*(pp.vortex^2)*[0 (1./gr.x(2:end)).^2] - 1 + (1/2)*pp.beta*(yb.^2))));
benergy = benergy + 2*pi*sum(gr.dx.*gr.x.*(-(1/2)*yb2prime));
benergy = benergy + 2*pi*sum(gr.dx.*(-(1/2)*ybprime));


end

% venergy = pi*sum(gr.dx.*gr.x.*(yvprime.^2 + [0 (yv(2:end)./gr.x(2:end)).^2] - yv.^2 + (1/2)*yv.^4));