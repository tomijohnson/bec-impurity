function  [Ei,ci]=solveimpurityeq(pp,np,gr,bb,yb,yi,na)

V = (1/2)*(1/pp.alpha)*(pp.omegat^2)*gr.x.^2+pp.etaab*pp.beta*(yb.^2)+pp.etaa*pp.beta*(na-1)*(yi.^2);

H = zeros(np.nb,np.nb);
for n = 1:np.nb
    for m = 1:np.nb
        H(n,m) = 2*pi*sum(gr.dx.*gr.x.*bb.J{n}.*bb.J{m}.*V);
    end
end

H = H + (1/2)*pp.alpha*diag(bb.z.^2)/np.maxx^2;

H = (H+H')/2;

[c,Ei] = eigs(H,1,'sa');

ci = c/sqrt(c'*c);

end