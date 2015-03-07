close all;
clear all;

nbasis = 20;

z = zerobess('J',0,nbasis);

ngrid = 5000;
maxx = 10;
a = 1/ngrid;
b = log(maxx/a+1)/(ngrid-1);
x = a*(exp(b*(linspace(1,ngrid+1,ngrid+1)-1))-1);
dx = diff(x);
x = x(1:end-1);

figure(1);
hold on;
for n = 1:nbasis
    J{n} = besselj(0,x*z(n)/maxx);
    J{n} = J{n}/sqrt(2*pi*sum(dx.*x.*(J{n}.^2)));
    plot(x,J{n});
end

for n = 1:nbasis
    for m = 1:nbasis
        G(n,m) = 2*pi*sum(dx.*x.*J{n}.*J{m});
    end
end

V = x.^2/2;

for n = 1:nbasis
    for m = 1:nbasis
        H(n,m) = 2*pi*sum(dx.*x.*J{n}.*J{m}.*V);
    end
end

H = H + (1/2)*diag(z.^2)/maxx^2;

H = (H+H')/2;

% opts.issym = 1;
% opts.isreal = 1;
% % opts.maxit = 1e3;
% % opts.p = 1e2;
% opts.v0 = yg1';

[c,E] = eigs(H,1,'sa');

c = c/sqrt(c'*c);

y = 0;
for n = 1:nbasis
    y = y + c(n)*J{n};
end

y = y*sign(y(1));

yg = ((1/pi)^(1/2))*exp(-x.^2/2);


figure(2);
plot(x,y);
hold on;
plot(x,yg,'r');