function Evg = getvenergytp2(pp,w)
	Evg = (pi/(2*pp.beta))*((1/2)*w^2-log(w^2)-(1-log(2)));
end