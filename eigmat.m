function T=eigmat(M,x,varargin)
    %parse for omega
    p=inputParser;
    p.addOptional('omega',1.0);
    
    parse(p,varargin{:});
    omega=p.Results.omega;
    
    if M>371
        error('double-precision quadrature does not converge for M>371');
    end
    
    psi0=exp(-(sqrt(omega)*x).^2/2)*(omega/pi)^(1/4);
    psi1=sqrt(2)*exp(-(sqrt(omega)*x).^2/2).*(sqrt(omega)*x)*(omega/pi)^(1/4);
    n=0:M-1;
    
    T=zeros(size(x.'*ones(1,M)));
    T(:,1)=psi0;
    T(:,2)=psi1;
    
    for m=1:M-2
        %size(sqrt(2/(n(m+2)))*(sqrt(omega)*x).*T(:,m+1)-sqrt(n(m+1)/n(m+2))*T(:,m))
        %size(T(:,m+2))
        T(:,m+2)=sqrt(2/(n(m+2)))*(sqrt(omega)*x.').*T(:,m+1)-sqrt(n(m+1)/n(m+2))*T(:,m);
    end
    
end