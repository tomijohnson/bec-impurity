clear all;
close all;



%Solving radial Schrodinger equation for harmonic oscillator in cylindrical coordinates

ngrid = 10000;
maxx = 10;
a = 1;
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

y1 = y1';




figure(1);
plot(x,y1);
hold on;
plot(x,yg1,'r');

y2 = y1.*(a*b*x.*exp(b*(linspace(1,ngrid+1,ngrid+1)-1))).^(-1/2);

figure(2);
plot(x,y2);
hold on;
plot(x,yg2,'r');

test = (H*(yg1'))./(yg1');



%linear grid
x = (linspace(1,ngrid+1,ngrid+1)-1)*10/(ngrid-1);
h = x(2)-x(1);


H = sparse(zeros(ngrid+1,ngrid+1));
    
%term containing potential
H = H + spdiags([fakepotential (1/2)*((x(2:end-1).^2)+(s^2-(1/4))./(x(2:end-1).^2)) fakepotential]',0,ngrid+1,ngrid+1);  

H = H + (1/(h)^2)*spdiags(ones(ngrid+1,ngrid+1)',0,ngrid+1,ngrid+1);

H = H - (1/2)*(1/(h)^2)*spdiags(ones(ngrid+1,ngrid+1)',1,ngrid+1,ngrid+1);

H = H - (1/2)*(1/(h)^2)*spdiags(ones(ngrid+1,ngrid+1)',-1,ngrid+1,ngrid+1);

H = (H+H')/2;

opts.issym = 1;
opts.isreal = 1;
opts.maxit = 5*1e3;
opts.p = 5*1e2;
opts.v0 = yg1';

[y1,Eg] = eigs(H,1,'sa');
    
% yg = (yg').*(x.^(-1/2));

y1 = y1';


figure(1);
plot(x,y1,'g');

y2 = y1.*(a*b*x.*exp(b*(linspace(1,ngrid+1,ngrid+1)-1))).^(-1/2);

figure(2);
plot(x,y2,'g');


test2 = (H*(yg1'))./(yg1');



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

opts.issym = 1;
opts.isreal = 1;
opts.maxit = 5*1e3;
opts.p = 5*1e2;
opts.v0 = yg1';

[y1,Eg] = eigs(H,1,'sa');
    
% yg = (yg').*(x.^(-1/2));

y1 = y1';


figure(1);
plot(x,y1,'g');

% y2 = y1.*(a*b*x.*exp(b*(linspace(1,ngrid+1,ngrid+1)-1))).^(-1/2);
% 
% figure(2);
% plot(x,y2,'g');


test2 = (H*(yg1'))./(yg1');

