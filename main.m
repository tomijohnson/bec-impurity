function main(job)
    addpath('gaussian','impurity','main','misc','tunnelling','twoparams1','twoparams2','bec');

    %----------------------------
    %    Load job details
    %----------------------------
    
    filename = 'jobdets.mat';
    load(filename,'pp','np','cp');
    display(['Loaded job details file - ' num2str(filename)]);
    job = str2num(job);
    display(['This is job ' num2str(job)]);
    display('The parameters are:');
    
    %----------------------------
    %    Physical parameters
    %----------------------------
    
    pp.alpha = pp.alpha(job);
    pp.beta = pp.beta(job);
    pp.etaa = pp.etaa(job);
    pp.etaab = pp.etaab(job);
    pp.omegat = pp.omegat(job);
    pp.namax = pp.namax(job);
    pp.vortex = pp.vortex(job);
    pp.becrad = pp.becrad(job);
    pp.dvals = pp.dvals(:,job);
    pp
    
    %----------------------------
    %    Numerical parameters
    %----------------------------
    
    np.maxx = np.maxx(job); %maximum radius we consider
    np.maxit1 = np.maxit(1,job);
    np.maxit2 = np.maxit(2,job);
    np.maxit3 = np.maxit(3,job);
    np=rmfield(np,{'maxit'});
    np.gridsize = np.gridsize(job); %number of points we wish to divide grid into
    np.acc = np.acc(job);
    np.nb = np.nb(job);
    np.nbr = np.nbr(job);
    np
    
    %----------------------------
    %    Computational parameters
    %----------------------------

    cp.dobg = cp.dobg(job); 
    cp.showbg = cp.showbg(job);
    cp.dobif = cp.dobif(job); 
    cp.showbif = cp.showbif(job);
    cp.dobi = cp.dobi(job);
    cp.showbi = cp.showbi(job);
    cp.doplot = cp.doplot(job);
%     cp.showtun = cp.showtun(job);
%     cp.dotun = cp.dotun(job);
    cp.savemin = cp.savemin(job);
    cp
    
    cp.filename = ['-' num2str(pp.alpha) '-' num2str(pp.beta) '-'  num2str(pp.etaa) '-' num2str(pp.etaab) '-' num2str(pp.omegat) '-' num2str(pp.namax) '-' num2str(pp.vortex) '-' num2str(pp.becrad) '-' num2str(np.maxx) '-' num2str(np.gridsize) '.mat'];
    
    %----------------------------
    %    Set up grid
    %----------------------------
    
    %nonuniform grid for emphasis around origin

    agrid = 1/np.gridsize;
    bgrid = log(np.maxx/agrid+1)/(np.gridsize-1);
    x = agrid*(exp(bgrid*(linspace(1,np.gridsize+1,np.gridsize+1)-1))-1);
    gr.dx = diff(x);
    gr.x = x(1:end-1);
    
    if (cp.dobg==1)||(cp.showbg==1)||(cp.dobif==1)||(cp.showbif==1)||(cp.dobi==1)||(cp.showbi==1)
        %----------------------------
        %    Solve with impurity within gaussian impurity approximation
        %----------------------------
        re = [];
        re = becgaussian(pp,np,cp,gr,re);
    end
    
    if (cp.dobif==1)||(cp.showbif==1)||(cp.dobi==1)||(cp.showbi==1)
        %----------------------------
        %    Solve with impurity without gaussian approx but with bec fixed
        %----------------------------

        re = becimpurityf(pp,np,cp,gr,re);
    end
    
    if (cp.dobi==1)||(cp.showbi==1)
        %----------------------------
        %    Solve with impurity without gaussian approx and without fixing
        %    bec
        %----------------------------

        re = becimpurity(pp,np,cp,gr,re);
    end
    
%     if (cp.dotun==1)||(cp.showtun==1)
%         %----------------------------
%         %    Caluculate tunelling of impurity between two wells
%         %----------------------------
%         re = tunnelling(pp,np,cp,gr,re);
%     end

    if (cp.doplot==1)
        doplot(pp,cp,gr,re);
    end

end