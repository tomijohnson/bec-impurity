addpath('../gaussian','../impurity','../main','../misc','../tunnelling','../twoparams1','../twoparams2','../vortex');

nmax = 3;

alpha = linspace(0.1,10,9);
beta = linspace(0.1,10,9);
etaa = linspace(0.1,10,9);
etaab = linspace(0.1,10,9);

% alpha = 1;
% beta = 1;
% etaa = 1;
% etaab = 10;

a = size(alpha,2);
b = size(beta,2);
c = size(etaa,2);
d = size(etaab,2);

alpha = repmat(alpha,[1 b*c*d]);
beta = repmat(beta,[a 1]);
beta = reshape(beta,[1 a*b]);
beta = repmat(beta,[1 c*d]);
etaa = repmat(etaa,[a*b 1]);
etaa = reshape(etaa,[1 a*b*c]);
etaa = repmat(etaa,[1 d]);
etaab = repmat(etaab,[a*b*c 1]);
etaab = reshape(etaab,[1 a*b*c*d]);

%     figure(1);
%     hold on;
%     plot(alpha);
%     plot(beta,'r');
%     plot(etaa,'g');
%     plot(etaab,'c');

sigma0 = 1;
w0 = 1;
rinit = [sigma0 w0];

successes1 = [];
successes2 = [];

for i = 1:size(alpha,2)
    
    pp.alpha = alpha(i);
    pp.beta = beta(i);
    pp.etaa = etaa(i);
    pp.etaab = etaab(i);
    
    valid1 = 0;
    valid2 = 0;
    
    for na = 1:nmax
        rmin1 = gettwoparams1(pp,na,rinit);
        Emin1(na) = getgenergy(pp,na,rmin1(1))+getvenergytp1(pp,rmin1(2))+getvgenergytp1(pp,na,rmin1(1),rmin1(2));
        rmin2 = gettwoparams2(pp,na,rinit);
        Emin2(na) = getgenergy(pp,na,rmin2(1))+getvenergytp2(pp,rmin2(2))+getvgenergytp2(pp,na,rmin2(1),rmin2(2));


        if na == 2
            if (Emin1(2) < 2*Emin1(1))
                valid1 = 1;
            end
            if (Emin2(2) < 2*Emin2(1))
                valid2 = 1;
            end
        end
        if (na > 2)&&(valid1 == 1)&&(valid1 ~= 2)&&(sum(na*Emin1(na) < (1:na-1).*Emin1(1:na-1),2) == 1)
            valid1 = 2;
            successes1 = [successes1 i];
        end
        if (na > 2)&&(valid2 == 1)&&(valid2 ~= 2)&&(sum(na*Emin2(na) < (1:na-1).*Emin2(1:na-1),2) == 1)
            valid2 = 2;
            successes2 = [successes2 i];
        end

        
    end
    
%     [Emin1' Emin2']'
end

successes1
successes2


%         if na == 1
%             re.U1(na) = 0;
%             re.U2(na) = 0;
%         else
%             re.U1(na) = (Emin1(na)-Emin1(1)*na)/(na*(na-1));
%             re.U2(na) = (Emin2(na)-Emin2(1)*na)/(na*(na-1));
%             if na == 2
%                 if (Emin1(2) < 2*Emin1(1))
%                     valid1 = 1;
%                 end
%                 if (Emin2(2) < 2*Emin2(1))
%                     valid2 = 1;
%                 end
%             end
%             if (na > 2)&&(valid1 == 1)&&(valid1 ~= 2)&&(sum(na*Emin1(na) < (1:na-1).*Emin1(1:na-1),2) == 1)
%                 valid1 = 2;
%                 successes1 = [successes1 i];
%             end
%             if (na > 2)&&(valid2 == 1)&&(valid2 ~= 2)&&(sum(na*Emin2(na) < (1:na-1).*Emin2(1:na-1),2) == 1)
%                 valid2 = 2;
%                 successes2 = [successes2 i];
%             end
%         end

% display('Showing successes');
% 
% successvalue = [];
% for i = successes
%     
%     a = alpha(i)/2;
%     b = etaab(i);
%     d = pi/(4*beta(i));
%     r = etaa(i)*beta(i)/(2*pi*alpha(i));
% 
%     param = [a b d r];
% 
%     [E_min,sol] = energyvsN(param,sol0);
%     dE1 = E_min(1)-E_min(2);
%     dE2 = E_min(3)-E_min(2);
%     valid = (dE1 > 0)&&(dE2 > 0);
% 
%     successvalue = [successvalue (dE1+dE2)/E_min(1)];
% end
% 
% [successvalue,ix] = sort(successvalue);
% successes = successes(ix);
% 
% for i = successes
% 
%     a = alpha(i)/2;
%     b = etaab(i);
%     d = pi/(4*beta(i));
%     r = etaa(i)*beta(i)/(2*pi*alpha(i));
% 
%     param = [a b d r];
% 
%     [E_min,sol] = energyvsN(param,sol0);
%     dE1 = E_min(1)-E_min(2);
%     dE2 = E_min(3)-E_min(2);
%     valid = (dE1 > 0)&&(dE2 > 0);
% 
%     if (valid==1)
%         param
%         [alpha(i) beta(i) etaa(i) etaab(i)]
%         sol
%         E_min
%     end
% end
