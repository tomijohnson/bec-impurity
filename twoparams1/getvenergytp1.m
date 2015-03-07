function Evg = getvenergytp1(pp,w)
	Evg = (pi/(2*pp.beta))*(w^2-log(w^2)-1);
end