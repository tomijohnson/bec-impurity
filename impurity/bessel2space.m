function y = bessel2space(np,bb,r,c)
y = 0;
for n = 1:np.nb
    y = y + c(n)*besselj(0,r*bb.z(n)/np.maxx)/bb.Jnorm(n);
end
y = y*sign(y(round(size(y,1)/2),round(size(y,2)/2)));
end