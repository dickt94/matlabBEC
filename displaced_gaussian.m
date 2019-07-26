function c0=displaced_gaussian(N,x0,width,nmodes,npaths)
    [x,w,T]=nfieldtrans(nmodes,2);
    invT=T';
    phi0=1/(pi^(1/4)*width^(1/2))*exp(-((x-x0)/width).^2/2);
    sqn0=invT*(w.*phi0);
    n0=N*sqn0.^2;
    phi0=repmat(2*pi*rand(1,npaths),nmodes,1);
    c0=zeros([nmodes+1 npaths],'gpuArray');
    c0(1:end-1,:)=sqrt(n0+1/2).*exp(1i*phi0);

end