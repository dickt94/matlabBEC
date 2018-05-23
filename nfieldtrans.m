%implements a subset of the functionality in
%https://github.com/AshtonSBradley/ProjectedGPE.jl/blob/master/src/nfieldtrans.jl,
%namely the Hermite basis transforms.

function [x,w,T] = nfieldtrans(M,n,varargin)
    %parse keyword args
    p=inputParser;
    p.addOptional('omega',1.0);
    
    parse(p,varargin{:});
    omega=p.Results.omega;
    
    if mod(n*M,2)==0
        K=n*M/2;
    else
        K=(n*M+1)/2;
    end
    [x,w]=gausshermite(K);
    w=exp(log(w)+x.^2)/sqrt(n*omega/2);
    T=eigmat(M,x/sqrt(n*omega/2),'omega',omega);
    
    x=x.';
    w=w.';
    
end

%memesss
