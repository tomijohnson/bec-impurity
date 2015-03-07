function lapy = laplacian(gr,y)

yex = [0 zeros(1,size(y,2)) 0; zeros(size(y,1),1) y zeros(size(y,1),1); 0 zeros(1,size(y,2)) 0];
dxex = [gr.dxx(1) gr.dxx gr.dxx(end)];
dyex = [gr.dyy(1) gr.dyy gr.dyy(end)];

lapy = 0;          

lapy = lapy + circshift(yex,[0 -1]).*(ones(gr.gridsizey+2,1)*( (1./dxex).*(1./((dxex + circshift(dxex,[0 -1]))/2))  ));
lapy = lapy + circshift(yex,[0 1]).*(ones(gr.gridsizey+2,1)*( (1./dxex).*(1./((dxex + circshift(dxex,[0 1]))/2))  ));
lapy = lapy - circshift(yex,[0 0]).*(ones(gr.gridsizey+2,1)*( (1./dxex).*(1./((dxex + circshift(dxex,[0 -1]))/2) + 1./((dxex + circshift(dxex,[0 1]))/2))  ));

lapy = lapy + circshift(yex,[-1 0]).*(( (1./dyex).*(1./((dyex + circshift(dyex,[0 -1]))/2))  )'*ones(1,gr.gridsizex+2));
lapy = lapy + circshift(yex,[1 0]).*(( (1./dyex).*(1./((dyex + circshift(dyex,[0 1]))/2))  )'*ones(1,gr.gridsizex+2));
lapy = lapy - circshift(yex,[0 0]).*(( (1./dyex).*(1./((dyex + circshift(dyex,[0 -1]))/2) + 1./((dyex + circshift(dyex,[0 1]))/2))  )'*ones(1,gr.gridsizex+2));

lapy = lapy(2:end-1,2:end-1);

end