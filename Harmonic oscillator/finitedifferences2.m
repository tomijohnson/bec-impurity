clear all;
close all;



%Solving radial Schrodinger equation for harmonic oscillator in cylindrical coordinates

ngrid = 10000;
maxx = 10;
a = 1;
b = log(maxx/a+1)/(ngrid-1);
x = a*(exp(b*(linspace(1,ngrid+1,ngrid+1)-1))-1);

fakepotential = 1e20;

s=0;


%linear symmetric grid
x = (linspace(1,ngrid+1,ngrid+1)-1)*10/(ngrid-1);
h = x(2)-x(1);
x = [-x x(2:end)];
ngrid = 2*ngrid+1;


H = sparse(zeros(ngrid,ngrid));
    
%term containing potential
H = H + spdiags((1/2)*((x.^2)+(s^2-(1/4))./(x.^2))',0,ngrid,ngrid);

H((ngrid-1)/2+1,(ngrid-1)/2+1) = fakepotential;

H = H + (1/(h)^2)*spdiags(ones(ngrid,ngrid)',0,ngrid,ngrid);

H = H - (1/2)*(1/(h)^2)*spdiags(ones(ngrid,ngrid)',1,ngrid,ngrid);

H = H - (1/2)*(1/(h)^2)*spdiags(ones(ngrid,ngrid)',-1,ngrid,ngrid);

H = (H+H')/2;

yg2 = ((1/pi)^(1/2))*exp(-x.^2/2);
yg1 = yg2.*(abs(x)).^(1/2);

opts.issym = 1;
opts.isreal = 1;
opts.maxit = 5*1e3;
opts.p = 5*1e2;
opts.v0 = yg1';

[y1,Eg] = eigs(H,1,'sa');

y1 = y1';


figure(1);
plot(x,y1,'g');

%known solution





test2 = (H*(yg1'))./(yg1');

