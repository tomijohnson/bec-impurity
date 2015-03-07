
alpha = linspace(0.1,10,8);
beta = linspace(0.1,10,8);
etaa = linspace(0.1,10,8);
etaab = linspace(0.1,10,8);

%     alpha = 3;
%     beta = 0.1;
%     etaa = 5;
%     etaab = 20;

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
sol0 = [sigma0 w0];

successes = [];

for i = 1:size(alpha,2)

    a = alpha(i)/2;
    b = etaab(i);
    d = pi/(4*beta(i));
    r = etaa(i)*beta(i)/(2*pi*alpha(i));

    param = [a b d r];

    [E_min,sol] = energyvsN(param,sol0);
    dE1 = E_min(1)-E_min(2);
    dE2 = E_min(3)-E_min(2);
    valid = (dE1 > 0)&&(dE2 > 0);

    if (valid==1)
        successes = [successes i];
    end
end

display('Showing successes');

successvalue = [];
for i = successes

    a = alpha(i)/2;
    b = etaab(i);
    d = pi/(4*beta(i));
    r = etaa(i)*beta(i)/(2*pi*alpha(i));

    param = [a b d r];

    [E_min,sol] = energyvsN(param,sol0);
    dE1 = E_min(1)-E_min(2);
    dE2 = E_min(3)-E_min(2);
    valid = (dE1 > 0)&&(dE2 > 0);

    successvalue = [successvalue (dE1+dE2)/E_min(1)];
end

[successvalue,ix] = sort(successvalue);
successes = successes(ix);

for i = successes

    a = alpha(i)/2;
    b = etaab(i);
    d = pi/(4*beta(i));
    r = etaa(i)*beta(i)/(2*pi*alpha(i));

    param = [a b d r];

    [E_min,sol] = energyvsN(param,sol0);
    dE1 = E_min(1)-E_min(2);
    dE2 = E_min(3)-E_min(2);
    valid = (dE1 > 0)&&(dE2 > 0);

    if (valid==1)
        param
        [alpha(i) beta(i) etaa(i) etaab(i)]
        sol
        E_min
    end
end
