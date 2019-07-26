function c0=displaced_gaussianhf(N,x0,width,nmodes)
    [x,w,T]=nfieldtrans(nmodes,2);
    invT=T';
    phi0=sqrt(N)/(pi^(1/4)*width^(1/2))*exp(-((x-x0)/width).^2/2);
    
    c0=invT*(w.*phi0);

end