function venergy = getbenergytest(pp,gr,yb)
%extract wavefunction and its derivative
yvext = [yb yb(end)];
yvprime = diff(yvext)./gr.dx;

%calculate energy
venergy = 2*pi*sum(gr.dx.*gr.x.*((1/2)*abs(yvprime).^2 + (yb.^2).*((1/2)*(pp.vortex^2)*[0 (1./gr.x(2:end)).^2] - 1 + (1/2)*pp.beta*(yb.^2))));
end