function create_job_details()

num_jobs = 1;
alljobs = ones(1,num_jobs);

%----------------------------
%    Physical parameters
%----------------------------

pp.alpha = 1*alljobs; % dimensionless parameters
pp.beta = 1*alljobs;
pp.etaa = 1*alljobs;
pp.etaab = 8*alljobs;
pp.omegat = 0*alljobs;

pp.namax = 3*alljobs; % maximum number of impurities at one site

pp.vortex = 1*alljobs; % is there a vortex or not?

pp.becrad = sqrt(1/pi)*alljobs; % where to we apply the boundary conditions for the BEC

pp.dvals = linspace(3,8,11)'*alljobs; % hopping parameters to consider

%----------------------------
%    Numerical parameters
%----------------------------

np.maxx = 50*alljobs; %maximum radius we consider
np.maxit = [20; 20; 20]*alljobs;
np.gridsize = 4998*alljobs; %number of points we wish to divide grid into
np.acc = 0.001*alljobs;
np.nb = 300*alljobs;
np.nbr = 150*alljobs;

%----------------------------
%    Computational parameters
%----------------------------

cp.dobg = 0*alljobs; 
cp.showbg = 1*alljobs;
cp.dobif = 0*alljobs; 
cp.showbif = 1*alljobs;
cp.dobi = 0*alljobs;
cp.showbi = 1*alljobs;
% cp.dotun = 0*alljobs;
% cp.showtun = 0*alljobs;
cp.doplot = 1*alljobs;
cp.savemin = 0*alljobs;

% % % ################################## EXAMPLE OF 2 #######################################
% a = size(chi,2);
% b = size(M,2);
% 
% chi = repmat(chi,[1 b]);
% M = reshape(repmat(M,[a 1]),[1 a*b]);
% 
% figure(1);
% hold off
% plot(chi/mean(chi));
% hold on;
% plot(M/mean(M),'r');
% 
% chi
% time.delt
% size(chi,2)
% num_jobs

% % #########################################################################

% ################################ EXAMPLE OF 3 #########################################
% a = size(pp.beta,2);
% b = size(pp.etaab,2);
% c = size(pp.etaa,2);
% 
% pp.beta = repmat(pp.beta,[1 b*c]);
% pp.etaab = reshape(repmat(pp.etaab,[a c]),[1 a*b*c]);
% pp.etaa = reshape(repmat(pp.etaa,[a*b 1]),[1 a*b*c]);
% 
% figure(1);
% hold off
% plot(pp.beta/mean(pp.beta));
% hold on;
% plot(pp.etaab/mean(pp.etaab),'r');
% plot(pp.etaa/mean(pp.etaa),'g');
% 
% size(pp.beta,2)
% num_jobs

% #########################################################################

filename = 'jobdets.mat';
display(['Saved job details file - ' num2str(filename)])
save(filename,'pp','np','cp');

end