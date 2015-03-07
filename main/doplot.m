function doplot(pp,cp,gr,re)
    for na=1:pp.namax
        figure(na);
        hold off; 
        if (cp.dobg==1)||(cp.showbg==1)||(cp.dobif==1)||(cp.showbif==1)||(cp.dobi==1)||(cp.showbi==1)
            dg = getgdensity(gr.x,re.bg.a(na));
            yg = dg.^(1/2);
            scale = max(yg);
            plot(gr.x,yg,'g');
            hold on;
            plot(gr.x,re.bg.yb{na}*scale/max(re.bg.yb{na}),'g-');
        end

        if (cp.dobif==1)||(cp.showbif==1)||(cp.dobi==1)||(cp.showbi==1)
            plot(gr.x,re.bif.yi{na},'r');
        end

        if (cp.dobi==1)||(cp.showbi==1)
            plot(gr.x,re.bi.yi{na});
            plot(gr.x,re.bi.yb{na}*scale/max(re.bi.yb{na}),'-');
        end
    end
%     if (cp.dotun==1)||(cp.showtun==1)
%         count = 0;
%         for na1=1:pp.namax
%             for na2=1:pp.namax
%                 count = count+1;
%                 dcount = 0;
%                 for d = pp.dvals'
%                     dcount = dcount+1;
%                     H(:,:,dcount) = re.H{dcount}{na1,na2};
%                     G(:,:,dcount) = re.G{dcount}{na1,na2};
%                     Htilde(:,:,dcount) = re.Htilde{dcount}{na1,na2};
%                 end
%                 
%                 t12 = -permute(Htilde(1,2,:),[3 1 2]);
%                 t21 = -permute(Htilde(2,1,:),[3 1 2]);
%                 
%                 
% 
%                 [estimates0,model0] = fitcurve0(pp.dvals,t12,[1 0.5]);
%                 [estimates,model] = fitcurve(pp.dvals,t12,[1 0.5]);
%                 
%                 [sse, FittedCurve] = model(estimates0);
%                 FittedCurve = FittedCurve*t12(end)/FittedCurve(end);
%                 figure(100+count);
%                 hold off;
%                 semilogy(pp.dvals,t12,'b');
%                 hold on;
%                 semilogy(pp.dvals,t21,'r');
%                 semilogy(pp.dvals,FittedCurve,'g');
%                 
% %                 figure(200+count);
% %                 hold off;
% %                 plot(pp.dvals,permute(Htilde(1,1,:),[3 1 2]),'b');
% %                 hold on;
% %                 plot(pp.dvals,permute(Htilde(2,2,:),[3 1 2]),'r');
% %                 plot(pp.dvals,(re.ienergy(na1)+re.vienergy(na1))*ones(size(permute(Htilde(2,2,:),[3 1 2]))),'g');
% %                 plot(pp.dvals,(re.ienergy(na2)+re.vienergy(na2))*ones(size(permute(Htilde(2,2,:),[3 1 2]))),'c');
% %                 
% 
%             end
%         end
%     end
end