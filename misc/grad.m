function gradysq = grad(gr,y)

yex = [0 zeros(1,size(y,2)) 0; zeros(size(y,1),1) y zeros(size(y,1),1); 0 zeros(1,size(y,2)) 0];
dxex = [gr.dxx(1) gr.dxx gr.dxx(end)];
dyex = [gr.dyy(1) gr.dyy gr.dyy(end)];

gradysq = 0;
gradysq = gradysq + ((circshift(yex,[0 -1])-circshift(yex,[0 1])).*(ones(gr.gridsizey+2,1)*(1./( dxex + circshift(dxex,[0 1]))))).^2;
gradysq = gradysq + ((circshift(yex,[-1 0])-circshift(yex,[1 0])).*((1./(dyex + circshift(dyex,[0 1])))'*(ones(1,gr.gridsizex+2)))).^2;
end