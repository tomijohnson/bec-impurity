function re = tunnelling(pp,np,cp,gr,re)

        display('--- Beginning calculation of tunnelling ---');
        
        %----------------------------
        %    Loop over d 
        %----------------------------
        dcount = 0;
        for d = pp.dvals'
            dcount = dcount + 1;
            filename = ['data/tunelling-' num2str(d) cp.filename];
            minfilename = ['data/mintunelling-' num2str(d) cp.filename];
            if ((exist(filename,'file') == 2)||(exist(['min' filename],'file') == 2))&&(cp.dotun~=1)
                if (exist(filename,'file') == 2)
                    load(filename,'tun');
                else
                    load(minfilename,'tun');
                end
                re.H{dcount} = tun.H;
                re.G{dcount} = tun.G;
                re.Htilde{dcount} = tun.Htilde;
                display(['--- Loaded previous tunnelling results from ' filename ' ---']);
            else
                display(['Well separation = ' num2str(d)]);
                %----------------------------
                %    Loop over occupations
                %----------------------------
                for na1 = 1:pp.namax
                    for na2 = 1:pp.namax 

                        display(['na1 = ' num2str(na1) '; na2 = ' num2str(na2)]);

    %                     gest(na1,na2,dcount) = 2*(a1*a2/(a1^2+a2^2))*exp(-(d^2)/(2*(a1^2+a2^2)))

                        %----------------------------
                        %    Build new grid
                        %----------------------------
                        a1 = re.a(na1);
                        a2 = re.a(na2);

                        %ygrid: nonuniform grid from -d/2 to d/2
                        averagea = (1/2)*(re.a(na1)+re.a(na2));
                        maxyy = 10*averagea;
                        np.gridsizey = 300;
                        agridy = 10*1/np.gridsizey;
                        bgridy = log(maxyy/agridy+1)/(np.gridsizey-1);
                        yy = agridy*(exp(bgridy*(linspace(1,np.gridsizey+1,np.gridsizey+1)-1))-1);
                        yy = [-fliplr(yy(2:end)) 0 yy(2:end)];
                        gr.dyy = diff(yy);
                        gr.dyy = (gr.dyy(1:end-1)+gr.dyy(2:end))/2;
                        gr.yy = yy(2:end-1);
                        gr.gridsizey = size(gr.yy,2);

                        %xgrid: nonuniform grid from -d to d, uniform in
                        %middle d and non-uniform either side
                        maxxx = 8*averagea;
                        np.gridsizex = 100;
                        agridx = 10*maxxx/np.gridsizex;
                        bgridx = log(maxxx/agridx+1)/(np.gridsizex-1);
                        xx = agridx*(exp(bgridx*(linspace(1,np.gridsizex+1,np.gridsizex+1)-1))-1);
                        xx = xx + d/2;
                        xx = [-fliplr(xx(2:end)) linspace(-d/2,d/2,d/(agridx*((exp(bgridx)-1)))+1) xx(2:end)];
                        gr.dxx = diff(xx);
                        gr.dxx = (gr.dxx(1:end-1)+gr.dxx(2:end))/2;
                        gr.xx = xx(2:end-1);
                        gr.gridsizex = size(gr.xx,2);

                        gr.xxx = ones(gr.gridsizey,1)*gr.xx;
                        gr.yyy = gr.yy'*ones(1,gr.gridsizex);

                        gr.dxdy = gr.dyy'*gr.dxx;

                        %----------------------------
                        %    Build merged vortex density on new grid
                        %----------------------------

                        %the precise way in which it is merged shouldn't matter
                        %too much

                        solv1 = re.solv{na1};
                        solv2 = re.solv{na2}; 

                        r1 = sqrt((gr.xxx+d/2).^2 +gr.yyy.^2);
                        r2 = sqrt((gr.xxx-d/2).^2 +gr.yyy.^2);

                        zerosind1 = find(r1 == 0);
                        zerosind2 = find(r2 == 0);

                        r1(zerosind1) = 1;
                        r2(zerosind2) = 1;

                        yv1 = deval(solv1,reshape(r1,[1 size(r1,1)*size(r1,2)]));
                        yv1 = yv1(1,:);
                        yv1 = reshape(yv1,[gr.gridsizey gr.gridsizex]);
                        yv1(zerosind1) = 0;
                        dv1 = yv1.^2;
                        clear yv1;
                        clear r1;

                        yv2 = deval(solv2,reshape(r2,[1 size(r2,1)*size(r2,2)]));
                        yv2 = yv2(1,:);
                        yv2 = reshape(yv2,[gr.gridsizey gr.gridsizex]);
                        yv2(zerosind2) = 0;
                        dv2 = yv2.^2;
                        clear yv2;
                        clear r2;

                        alpha = 5;
                        dv = (dv1.*exp(-alpha*gr.xxx) + dv2.*exp(alpha*gr.xxx))./(2*cosh(alpha*gr.xxx));
                        clear dv1;
                        clear dv2;

                        %----------------------------
                        %    Build impurity wavefunctions on new grid
                        %----------------------------
                        r1 = sqrt((gr.xxx+d/2).^2 +gr.yyy.^2);
                        r2 = sqrt((gr.xxx-d/2).^2 +gr.yyy.^2);

    %                     yi1 = exp(-r1.^2/(2*a1^2))*(pi*a1^2)^(-1/2);
    %                     yi2 = exp(-r2.^2/(2*a2^2))*(pi*a2^2)^(-1/2);
    %                     sum(sum(dxdy.*yi1.^2))
    %                     sum(sum(dxdy.*yi2.^2))
    %                     g(na1,na2,dcount) = sum(sum(dxdy.*yi1.*yi2))

                        bb = buildbasis(np,gr);

                        yi1 = bessel2space(np,bb,r1,re.ci{na1});
                        yi2 = bessel2space(np,bb,r2,re.ci{na2});

                        clear r1;
                        clear r2;


                        %----------------------------
                        %    Extract interesting quantities
                        %----------------------------       

                        %norms (to check that grid is accurate enough)
                        n1 = sum(sum(gr.dxdy.*yi1.^2));
                        n2 = sum(sum(gr.dxdy.*yi2.^2));

                        %overlaps
                        g12 = sum(sum(gr.dxdy.*yi1.*yi2));
                        g21 = sum(sum(gr.dxdy.*yi2.*yi1));

                        %diagonals
                        e1 = sum(sum(gr.dxdy.*(-(1/2)*yi1.*laplacian(gr,yi1)+pp.etaab*pp.beta*yi1.*dv.*yi1)));
                        e2 = sum(sum(gr.dxdy.*(-(1/2)*yi2.*laplacian(gr,yi2)+pp.etaab*pp.beta*yi2.*dv.*yi2)));

                        %alternative diagonals using grad
    %                     e1alt = sum(sum(gr.dxdy.*((1/2)*gradsq(gr,yi1)+pp.etaab*pp.beta*dv.*yi1.^2)))
    %                     e2alt = sum(sum(gr.dxdy.*((1/2)*gradsq(gr,yi2)+pp.etaab*pp.beta*dv.*yi2.^2)))

                        %perhaps worth checking that not too far off what is
                        %expected without the other well there
    %                     e1altalt = re.ienergy(na1)+re.vienergy(na1)
    %                     e1altalt = re.ienergy(na1)+re.vienergy(na1)

                        %tunnelling as overlap with hamiltonian
                        t12 = sum(sum(gr.dxdy.*(-(1/2)*yi1.*laplacian(gr,yi2)+pp.etaab*pp.beta*yi1.*dv.*yi2)));
                        t21 = sum(sum(gr.dxdy.*(-(1/2)*yi2.*laplacian(gr,yi1)+pp.etaab*pp.beta*yi2.*dv.*yi1)));

                        %alternative tunnelling as difference between symmetric and
                        %antisymmetric energies (good approx if overlap/tunnelling small)
    %                     ys = (yi1+yi2)/sqrt(2);
    %                     ya = (yi1-yi2)/sqrt(2);
    %                     
    %                     Es = sum(sum(gr.dxdy.*(-(1/2)*ys.*laplacian(gr,ys)+pp.etaab*pp.beta*dv.*ys.^2)));
    %                     Ea = sum(sum(gr.dxdy.*(-(1/2)*ya.*laplacian(gr,ya)+pp.etaab*pp.beta*dv.*ya.^2)));
    %                     
    %                     Esalt = sum(sum(gr.dxdy.*((1/2)*gradsq(gr,ys)+pp.etaab*pp.beta*dv.*ys.^2)));
    %                     Eaalt = sum(sum(gr.dxdy.*((1/2)*gradsq(gr,ya)+pp.etaab*pp.beta*dv.*ya.^2)));
    %                     
    %                     t12alt = (Es-Ea)/2;
    %                     t12altalt = (Esalt-Eaalt)/2;

                        %construct Hamiltonian in non-orthogonal basis
                        H = [e1 t21; t12 e2];
                        G = [n1 g21; g12 n2];

                        %move to Hamiltonian in orthogonal basis
                        Htilde = (G^(-1/2))'*H*(G^(-1/2));

                        %construct antisymmetric/symmetric combinations of
                        %orthogonal basis
    %                     S = G^(-1/2);
    %                     ystilde = ((S(1,1)+S(1,2))*yi1+(S(2,1)+S(2,2))*yi2)/sqrt(2);
    %                     yatilde = ((S(1,1)-S(1,2))*yi1+(S(2,1)-S(2,2))*yi2)/sqrt(2);  
    %                     nstilde = sum(sum(gr.dxdy.*ystilde.^2))
    %                     natilde = sum(sum(gr.dxdy.*yatilde.^2))
    %                     Estildealt = sum(sum(gr.dxdy.*(-(1/2)*ystilde.*laplacian(gr,ystilde)+pp.etaab*pp.beta*dv.*ystilde.^2)))
    %                     Eatildealt = sum(sum(gr.dxdy.*(-(1/2)*yatilde.*laplacian(gr,yatilde)+pp.etaab*pp.beta*dv.*yatilde.^2)))
    %                     ttildealt = (Estildealt-Eatildealt)/2;

                        %----------------------------
                        %    Display results
                        %----------------------------

                        display(['Tunelling = ' num2str(Htilde(1,2)) ', ' num2str(Htilde(2,1)) '; Energies = ' num2str(Htilde(1,1)) ', ' num2str(Htilde(2,2))]);

                        %----------------------------
                        %    Store results
                        %----------------------------
                        tun.H{na1,na2} = H;
                        tun.G{na1,na2} = G;
                        tun.Htilde{na1,na2} = Htilde;
                        
                        re.H{dcount}{na1,na2} = H;
                        re.G{dcount}{na1,na2} = G;
                        re.Htilde{dcount}{na1,na2} = Htilde;
                        
                        %----------------------------
                        %    Save result
                        %----------------------------
                        if cp.savemin == 1;
                            re2 = re;
                            re=rmfield(re,{'yv','yvg','yvf','solv','solvg','solvf','yi','ci'});
                            save(minfilename,'pp','np','re','tun');
                            display(['Results saved to ' minfilename])
                            re = re2;
                        else
                            save(filename,'pp','np','re','tun');
                            display(['Results saved to ' filename])
                        end
                        

                        %----------------------------
                        %    Plot results
                        %----------------------------

    %                     figure(1);
    %                     surf(gr.xx,gr.yy,yi1);
    %                     shading interp;
    %                     
    %                     figure(2);
    %                     surf(gr.xx,gr.yy,yi2);
    %                     shading interp;
    %                     
    %                     figure(3);
    %                     surf(gr.xx,gr.yy,dv);
    %                     shading interp;

    %                     figure(4);
    %                     surf(gr.xx,gr.yy,ys);
    %                     shading interp;
    %                     
    %                     figure(5);
    %                     surf(gr.xx,gr.yy,ya);
    %                     shading interp;

    %                     figure(6);
    %                     surf(gr.xx,gr.yy,ystilde);
    %                     shading interp;
    %                     
    %                     figure(7);
    %                     surf(gr.xx,gr.yy,yatilde);
    %                     shading interp;
 
                    end
                end
            end
        end

%         %----------------------------
%         %    Display result
%         %----------------------------
%         display('Summary of occupation dependence of total and interaction energies');
%         display('n, E_n, U_n, a, h');    
%         display([num2str(0) ', ' num2str(0) ', ' num2str(0) ', na, ' num2str(re.densityhole0)]);
%         for na = 1:pp.namax
%             display([num2str(na) ', ' num2str(re.totenergyg(na)) ', ' num2str(re.Ug(na)) ', ' num2str(re.a(na)) ', ' num2str(re.densityhole(na))]);
%         end

        display('--- Completed calculation of tunnelling ---');

end