function bb = buildbasis(np,gr)
    bb.z = zerobess('J',0,np.nb);
    for n = 1:np.nb
        J{n} = besselj(0,gr.x*bb.z(n)/np.maxx);
        bb.Jnorm(n) = sqrt(2*pi*sum(gr.dx.*gr.x.*(abs(J{n}).^2)));
        bb.J{n} = J{n}/bb.Jnorm(n);
    end
end