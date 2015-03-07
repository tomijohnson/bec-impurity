function re = becimpurityf(pp,np,cp,gr,re)
    filename = ['data/bif' cp.filename];
    minfilename = ['data/minbif' cp.filename];
    
    if (exist(filename,'file') == 2)&&(cp.dobif~=1)
        load(filename,'re');
        display(['--- Loaded previous calculation with fixed bec from ' filename ' ---']);
    else
        display('--- Beginning impurity calculation with fixed bec ---');  
        %----------------------------
        %    Set up Bessel basis functions
        %----------------------------

        bb = buildbasis(np,gr);

        for na = 1:pp.namax 

            display(['na = ' num2str(na)])

            %----------------------------
            %    Solves for impurity using iterative procedure
            %    without gaussian approximation, but fixed BEC
            %----------------------------

            %use gaussian results as initial guesses
            dg = getgdensity(gr.x,re.bg.a(na));
            yg = dg.^(1/2);
            yi = yg;
            yb = re.bg.yb{na};

            %begin iterative procedure
            itcount = 0;
            iterate = 1;      
            while iterate == 1
                itcount = itcount+1;

                %----------------------------
                %    Solves for impurity with fixed BEC and impurity self-density
                %----------------------------

                [Eiit(itcount+1),ci]=solveimpurityeq(pp,np,gr,bb,yb,yi,na);
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
                %    Check for convergence
                %----------------------------
                if itcount > 2
                    iterate = checkconvergence(Eiit,np.acc);
                end

                if itcount > np.maxit1
                    iterate = 0;
                end
            end

            %----------------------------
            %    Extract interesting quantities
            %----------------------------

            Ei = Eiit(end);
            clear Eiit;

            %store energies
            re.bif.benergy(na) = getbenergy(pp,gr,yb)-re.Eb0;
            re.bif.ienergy(na) = getienergy(pp,gr,yi,na);
            re.bif.bienergy(na) = getbienergy(pp,gr,yb,yi,na);
            re.bif.totenergy(na) = re.bif.benergy(na)+re.bif.ienergy(na)+re.bif.bienergy(na);


            if na == 1
                re.bif.U(na) = 0;
            else
                re.bif.U(na) = (re.bif.totenergy(na)-re.bif.totenergy(1)*na)/(na*(na-1));
            end

            %store everything else
            re.bif.yi{na} = yi;
            re.bif.ci{na} = ci;

            %----------------------------
            %    Display result
            %----------------------------

            display(['Tot energy = ' num2str(re.bif.totenergy((na))) '; Int energy = ' num2str(re.bif.U(na))]);

        end

        %----------------------------
        %    Display result
        %----------------------------
        display('Summary of occupation dependence of total and interaction energies');
        display('n, E_n, U_n');    
        display([num2str(0) ', ' num2str(0) ', ' num2str(0)]);
        for na = 1:pp.namax
            display([num2str(na) ', ' num2str(re.bif.totenergy(na)) ', ' num2str(re.bif.U(na))]);
        end

        %----------------------------
        %    Save result
        %----------------------------
        if cp.savemin == 1;
            re2 = re;
            re.bif=rmfield(re.bif,{'yi','ci'});
            save(minfilename,'pp','np','re.bif');
            display(['Results saved to ' minfilename])
            re = re2;
        else
            save(filename,'pp','np','re');
            display(['Results saved to ' filename])
        end
        display('--- Completed impurity calculation with fixed bec ---');
    end
end