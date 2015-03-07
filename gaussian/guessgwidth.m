function a = guessgwidth(gr,pp,yb)
a = (2*(((yb(3)-yb(2))/(gr.x(3)-gr.x(2)))^2)*pp.etaab*pp.beta/pp.alpha)^(-1/4);
end