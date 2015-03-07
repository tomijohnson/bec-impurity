%script to test all of the functions used
close all;
clear all;

%----------------------------
%    Physical parameters
%----------------------------

pp.s = 1;
pp.alpha = 2.3;
pp.beta = 0.9;
pp.etaab = 2;
pp.namax = 1;
pp.etaa = 3;

%----------------------------
%    Numerical parameters
%----------------------------

np.maxx = 50; %maximum radius we consider
np.maxit1 = 20;
np.maxit2 = 20;
np.gridsize = 4998; %number of points we wish to divide grid into
np.acc = 0.001;
np.nb = 200;

%----------------------------
%    Set up grid
%----------------------------

%nonuniform grid for emphasis around origin

agrid = 1/np.gridsize;
bgrid = log(np.maxx/agrid+1)/(np.gridsize-1);
x = agrid*(exp(bgrid*(linspace(1,np.gridsize+1,np.gridsize+1)-1))-1);
gr.dx = diff(x);
gr.x = x(1:end-1);

%----------------------------
%    Test vortex related functions
%----------------------------

%run function with no impurity
yv=solvevortexeqg(pp,gr,0,1,0);

%plot to check it all looks sensible, going to 1/sqrt(beta) at the boundary
figure(1);
plot(gr.x,yv);
hold on;
plot(gr.x,(1/sqrt(pp.beta))*gr.x./(gr.x.^2+1).^(1/2),'r'); %plot against well known approximation function
hold off;

%test to see if laplacian replacement ok
venergy = getvenergy(pp,gr,yv)
venergy = getvenergytest(pp,gr,yv)

%see how accurate replacement of the chemical potential by ng is by comparing densities.
Nb = 2*pi*sum(gr.dx.*gr.x.*(yv.^2));
density = Nb/(pi*np.maxx^2)
1/pp.beta

%find a way to explicitly estimate the chemical potential

%check density hole maxes sense
vdensityhole = getvdensityhole(pp,gr,yv);
effvrad = sqrt(pp.beta*vdensityhole/pi)

%----------------------------
%    Test gaussian related functions
%----------------------------
omega = 1;

yv = sqrt(omega^2/(2*pp.alpha*pp.etaab*pp.beta))*gr.x;
a = guessgwidth(gr,pp,yv)
a = solvegaussianeq(pp,np,gr,yv,1,a)
a = sqrt(pp.alpha/omega)

%----------------------------
%    Test impurity related functions
%----------------------------

%buld bessel basis
bb = buildbasis(np,gr);

%check normalised
G = zeros(np.nb,np.nb);
for n = 1:np.nb
    for m = 1:np.nb
        G(n,m) = 2*pi*sum(gr.dx.*gr.x.*bb.J{n}.*bb.J{m});
    end
end
max(max(abs(G-eye(np.nb))))

%check energies
for n = 1:np.nb
        FreeEnergy(n) = getienergy(pp,gr,bb.J{n},1);
        FreeEnergy2(n) = (1/2)*pp.alpha*bb.z(n)^2/np.maxx^2;
end
[FreeEnergy; FreeEnergy2]

%give harmonic oscillator/gaussian test to check wavefunction and energies of full calculation
dg = getgdensity(gr.x,a);
yg = dg.^(1/2);
yi = yg;
[Ei,ci]=solveimpurityeq(pp,np,gr,bb,yv,yi,1);
yi = 0;
for n = 1:np.nb
    yi = yi + ci(n)*bb.J{n};
end

%extract solution
yi = yi*sign(yi(1));

%plot against known gaussian soltuion
figure(2);
plot(gr.x,yi);
hold on;
plot(gr.x,yg,'r');
hold off;

%test impurity energies
ienergy = getienergy(pp,gr,yi,1)
ienergy = getienergy(pp,gr,yg,1)
genergy = getgenergy(pp,1,a)
omega/2

%test interaction energies
vienergy = getvienergy(pp,gr,yv,yi,1)
vienergy = getvienergy(pp,gr,yv,yg,1) %probably differ slightly as yi is not exactly gaussian
vgenergy = getvgenergy(pp,gr,yv,1,a)   
omega/2

%test total energies
totenergy = ienergy+vienergy
totenergyg = genergy+vgenergy
omega
Ei

%test to see if laplacian replacement ok
ienergy = getienergy(pp,gr,yi,1)
ienergytest = getienergytest(pp,gr,yi,1)

%now go on to test numerically and by checking the analytics that the na>1
%case is ok. easy to test for etaab = 0, but here we want to see more
%generally

%test impurity energies
ienergy = getienergy(pp,gr,yi,2)
genergy = getgenergy(pp,2,a)

%test interaction energies
vienergy = getvienergy(pp,gr,yv,yi,2)
vgenergy = getvgenergy(pp,gr,yv,2,a)

%test total energies
totenergy = ienergy+vienergy
totenergyg = genergy+vgenergy

%in the three tests above, we have used y1 that is the n=1 solution. maybe
%now go onto calculat the n=2 solution and then perform the test again, and
%compare with Ei

itcount = 0;
iterate = 1;      
while iterate == 1
    itcount = itcount+1;

    %----------------------------
    %    Solves for impurity with vortex and fixed impurity self-density
    %----------------------------

    [Eiit(itcount+1),ci]=solveimpurityeq(pp,np,gr,bb,yv,yi,2);
    yi = 0;
    for n = 1:np.nb
        yi = yi + ci(n)*bb.J{n};
    end
    yi = yi*sign(yi(1));

    %----------------------------
    %    Display result
    %----------------------------

    display(['Iteration number = ' num2str(itcount) '; Impurity eigenvalue = ' num2str(Eiit(itcount+1))]);

    %----------------------------
    %    Check for convergence (maybe energy is not so good, as
    %    won't be sensitive to the tail)
    %----------------------------
    if itcount > 2
        iterate = checkconvergence(Eiit,np.acc);
    end

    if itcount > np.maxit1
        iterate = 0;
    end
end

%test impurity energies
ienergy = getienergy(pp,gr,yi,2)
genergy = getgenergy(pp,2,a) %should be larger, as a is too small a width

%test interaction energies
vienergy = getvienergy(pp,gr,yv,yi,2)
vgenergy = getvgenergy(pp,gr,yv,2,a) %should be smaller, as a is too small a width

%test total energies
totenergy = ienergy+vienergy
totenergyg = genergy+vgenergy %should be larger, as will have energy larger than ground state

%impurity plus interaction
ienergy
Ei = Eiit(end) %should be larger, as takes self-interaction energy into account twice


% w = 1;
% yv = (1/sqrt(pp.beta))*gr.x./(gr.x.^2+w^2).^(1/2);
getvtpenergy(pp,np,1)
getvtpenergy(pp,np,0.5)
getvtpenergy(pp,np,1.5)

getdeltvtpenergy(pp,0.5)
getvtpenergy(pp,np,0.5) - getvtpenergy(pp,np,1)

getdeltvtpenergy(pp,1.5)
getvtpenergy(pp,np,1.5) - getvtpenergy(pp,np,1)