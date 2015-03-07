clear all;
close all;



%Solving radial Schrodinger equation for harmonic oscillator in cylindrical coordinates

ngrid = 10000;
maxx = 10;
a = 0.00001;
b = log(maxx/a+1)/(ngrid-1);
x = a*(exp(b*(linspace(1,ngrid+1,ngrid+1)-1))-1);

fakepotential = 1e6;

s=0;

H = sparse(zeros(ngrid+1,ngrid+1));
    
%term containing potential
H = H + spdiags([fakepotential (1/2)*((x(2:end-1).^2)+(s^2-(1/4))./(x(2:end-1).^2)) fakepotential]',0,ngrid+1,ngrid+1);     

H = H + (1/(a*b)^2)*spdiags(1./exp(2*b*(linspace(1,ngrid+1,ngrid+1)-2))',0,ngrid+1,ngrid+1); 

H = H - ((1/2)*(exp(b)/cosh(b/2))*(1/(a*b)^2)*spdiags(1./exp(2*b*linspace(1,ngrid+1,ngrid+1))',-1,ngrid+1,ngrid+1))';

H = H - (1/2)*(exp(-b)/cosh(b/2))*(1/(a*b)^2)*spdiags(1./exp(2*b*(linspace(1,ngrid+1,ngrid+1)-1))',-1,ngrid+1,ngrid+1);

H = (H+H')/2;



%known solution

yg2 = ((1/pi)^(1/2))*exp(-x.^2/2);
yg1 = yg2.*(a*b*x.*exp(b*(linspace(1,ngrid+1,ngrid+1)-1))).^(1/2);

%find solution

opts.issym = 1;
opts.isreal = 1;
% opts.maxit = 1e3;
% opts.p = 1e2;
opts.v0 = yg1';

[y1,Eg] = eigs(H,1,'sa',opts);
    
% yg = (yg').*(x.^(-1/2));


y1 = y1'/sqrt(sum(y1'*y1));
y2 = y1.*(a*b*x.*exp(b*(linspace(1,ngrid+1,ngrid+1)-1))).^(-1/2);

dx = diff(x);
x = x(1:end-1);
yg1 = yg1(1:end-1);
yg2 = yg2(1:end-1);
y1 = y1(1:end-1);
y2 = y2(1:end-1);
y2(1) = y2(2);

y2 = y2/sqrt(2*pi*sum(dx.*x.*abs(y2).^2));


figure(1);
plot(x,y1);
hold on;
plot(x,yg1,'r');





figure(2);
plot(x,y2);
hold on;
plot(x,yg2,'r');

% test = (H*(yg1'))./(yg1');


% 2*pi*sum(dx.*x.*abs(yg2).^2)
