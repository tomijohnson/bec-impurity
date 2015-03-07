function re = becimpurity(pp,np,cp,gr,re)
    filename = ['data/bi' cp.filename];
    minfilename = ['data/minbi' cp.filename];
    
    if (exist(filename,'file') == 2)&&(cp.dobi~=1)
        load(filename,'re');
        display(['--- Loaded previous calculation from ' filename ' ---']);
    else
        display('--- Beginning full impurity calculation ---');  
        %----------------------------
        %    Set up Bessel basis functions
        %----------------------------

        bb = buildbasis(np,gr);
        bbr = buildbasis(np,gr);

        for na = 1:pp.namax 

            display(['na = ' num2str(na)])

            %----------------------------
            %    Solves for vortex with impurity using iterative procedure
            %    without gaussian approximation
            %----------------------------

            %use previous results as initial guesses
            yi = re.bif.yi{na};
            yb = re.bg.yb{na};

            %begin iterative procedure
            itcount = 0;
            iterate = 1;    
            Etotit(1) = re.bif.totenergy(na);
            Ei = re.bif.ienergy(na)+re.bif.bienergy(na); %this is not actually going to be equal to Ei, due to the factor of 2 in the interaction energy
            while iterate == 1
                itcount = itcount+1;


                %----------------------------
                %    Solves for impurity with vortex and fixed impurity self-density
                %----------------------------

                %begin iterative procedure
                itcount2 = 0;
                iterate2 = 1;      
                Eiit = Ei;
                while iterate2 == 1
                    itcount2 = itcount2+1;
                    [Eiit(itcount2+1),ci]=solveimpurityeq(pp,np,gr,bb,yb,yi,na);
                    yi = 0;
                    for n = 1:np.nb
                        yi = yi + ci(n)*bb.J{n};
                    end
                    yi = yi*sign(yi(1));

                    %----------------------------
                    %    Display result
                    %----------------------------

                    display(['    Imp only it no = ' num2str(itcount2) '; Impurity eigenvalue = ' num2str(Eiit(itcount2+1))]);


                    %----------------------------
                    %    Check for convergence (maybe energy is not so good, as
                    %    won't be sensitive to the tail)
                    %----------------------------
                    if itcount2 > 2
                        iterate2 = checkconvergence(Eiit,np.acc);
                    end

                    if itcount2 > np.maxit3
                        iterate2 = 0;
                    end
                end

                %----------------------------
                %    Solves for vortex with impurity
                %----------------------------

                [yb,solb]=solvebeceq(pp,np,gr,bbr,ci,na);

                %----------------------------
                %    Display result
                %----------------------------

                Ei = Eiit(end);
                Etotit(itcount+1) = getbenergy(pp,gr,yb)-re.Eb0+Ei;
                % Alternatively Etotit(itcount+1) = getbenergy(pp,gr,yb)-re.Eb0+getienergy(pp,gr,yi,na)+getbienergy(pp,gr,yb,yi
                clear Eiit;
                display(['Iteration number = ' num2str(itcount) '; Total energy = ' num2str(Etotit(itcount+1))]);

                %----------------------------
                %    Check for convergence
                %----------------------------
                if itcount > 2
                    iterate = checkconvergence(Etotit,np.acc);
                end

                if itcount > np.maxit1
                    iterate = 0;
                end
            end

            %----------------------------
            %    Extract interesting quantities
            %----------------------------

            %store energies
            re.bi.benergy(na) = getbenergy(pp,gr,yb)-re.Eb0;
            re.bi.ienergy(na) = getienergy(pp,gr,yi,na);
            re.bi.bienergy(na) = getbienergy(pp,gr,yb,yi,na);
            re.bi.totenergy(na) = re.bi.benergy(na)+re.bi.ienergy(na)+re.bi.bienergy(na);


            if na == 1
                re.bi.U(na) = 0;
            else
                re.bi.U(na) = (re.bi.totenergy(na)-re.bi.totenergy(1)*na)/(na*(na-1));
            end

            %store everything else
            re.bi.yi{na} = yi;
            re.bi.ci{na} = ci;
            re.bi.yb{na} = yb;

            %----------------------------
            %    Display result
            %----------------------------

            display(['Tot energy = ' num2str(re.bi.totenergy((na))) '; Int energy = ' num2str(re.bi.U(na))]);

        end

        %----------------------------
        %    Display result
        %----------------------------
        display('Summary of occupation dependence of total and interaction energies');
        display('n, E_n, U_n');    
        display([num2str(0) ', ' num2str(0) ', ' num2str(0)]);
        for na = 1:pp.namax
            display([num2str(na) ', ' num2str(re.bi.totenergy(na)) ', ' num2str(re.bi.U(na))]);
        end

        %----------------------------
        %    Save result
        %----------------------------
        if cp.savemin == 1;
            re2 = re;
            re.bi=rmfield(re.bi,{'yb','yi','ci'});
            save(minfilename,'pp','np','re.bi');
            display(['Results saved to ' minfilename])
            re = re2;
        else
            save(filename,'pp','np','re');
            display(['Results saved to ' filename])
        end
        display('--- Completed full impurity calculation ---');
        
    end
end