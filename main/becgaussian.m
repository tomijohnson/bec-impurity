function re = becgaussian(pp,np,cp,gr,re)
    filename = ['data/bg' cp.filename];
    minfilename = ['data/minbg' cp.filename];
    
    if (exist(filename,'file') == 2)&&(cp.dobg~=1)
        load(filename,'re');
        display(['--- Loaded previous calculation within gaussian approx from ' filename ' ---']);
    else
        display('--- Beginning calculation within gaussian approx ---');

        %----------------------------
        %    Solves for bec without impurity first to give location of zero
        %    energy and estimate width
        %----------------------------

        
        %this solves the bec equation without the impurity
        display(['na = ' num2str(0)])
        na = 0;
        [yb,solb]=solvebeceqg(pp,np,gr,0,1); %arbitrary value for a is given as does not contribute for na = 0
        
        %get zero energy
        re.Eb0 = getbenergy(pp,gr,yb);
        re.densityhole0 = getbdensityhole(pp,gr,yb);

        %make guess for impurity width
        a = 1;

        display(['BEC energy ' num2str(re.Eb0) ' set as zero; Impurity width guess = ' num2str(a) '; V hole unpeturbed = ' num2str(re.densityhole0)]);

        %----------------------------
        %    Solves for BEC with impurity using iterative procedure
        %----------------------------
        for na = 1:pp.namax 

            display(['na = ' num2str(na)])

            %pick an initial value for the impurity width, but stop it from
            %being too small
            ait(1) = max(a,0.2);

            %begin iterative procedure
            itcount = 0;
            iterate = 1;      
            while iterate == 1
                itcount = itcount+1;

                %----------------------------
                %    Solves for BEC with impurity
                %----------------------------

                [yb,solb]=solvebeceqg(pp,np,gr,na,ait(itcount));

                %----------------------------
                %    Solves for impurity with vortex
                %----------------------------

                ait(itcount+1)=solvegaussianeq(pp,np,gr,yb,na,ait(itcount));

                %----------------------------
                %    Display result
                %----------------------------

                display(['Iteration number = ' num2str(itcount) '; Impurity width = ' num2str(ait(itcount+1))]);

                %----------------------------
                %    Check for convergence
                %----------------------------
                if itcount > 2
                    iterate = checkconvergence(ait,np.acc);
                end

                if itcount > np.maxit1
                    iterate = 0;
                end
            end

            %----------------------------
            %    Extract interesting quantities
            %----------------------------

            a = ait(end);
            clear ait;

            %store energies
            re.bg.benergy(na) = getbenergy(pp,gr,yb)-re.Eb0;
            re.bg.ienergy(na) = getgenergy(pp,na,a);
            re.bg.bienergy(na) = getbgenergy(pp,gr,yb,na,a);
            re.bg.totenergy(na) = re.bg.benergy(na)+re.bg.ienergy(na)+re.bg.bienergy(na);

            if na == 1
                re.bg.U(na) = 0;
            else
                re.bg.U(na) = (re.bg.totenergy(na)-re.bg.totenergy(1)*na)/(na*(na-1));
            end

            %store everything else
            re.bg.a(na) = a;
            re.bg.yb{na} = yb;
            re.bg.densityhole(na) = getbdensityhole(pp,gr,yb);
            

            %----------------------------
            %    Display result
            %----------------------------

            display(['Tot energy = ' num2str(re.bg.totenergy(na)) '; Int energy = ' num2str(re.bg.U(na)) '; imp width = ' num2str(re.bg.a(na)) '; bec hole = ' num2str(re.bg.densityhole(na))]);

        end

        %----------------------------
        %    Display result
        %----------------------------
        display('Summary of occupation dependence of total and interaction energies');
        display('n, E_n, U_n, a, h');    
        display([num2str(0) ', ' num2str(0) ', ' num2str(0) ', na, ' num2str(re.densityhole0)]);
        for na = 1:pp.namax
            display([num2str(na) ', ' num2str(re.bg.totenergy(na)) ', ' num2str(re.bg.U(na)) ', ' num2str(re.bg.a(na)) ', ' num2str(re.bg.densityhole(na))]);
        end

        %----------------------------
        %    Save result
        %----------------------------
        if cp.savemin == 1;
            re2 = re;
            re.bg=rmfield(re.bg,{'yb'});
            save(minfilename,'pp','np','re.bg');
            display(['Results saved to ' minfilename])
            re = re2;
        else
            save(filename,'pp','np','re');
            display(['Results saved to ' filename])
        end
        display('--- Completed calculation within gaussian approx ---');
        
    end
end