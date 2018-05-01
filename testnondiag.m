function [S,T]=testnondiag()
%implements the first non-diagonal noise example given at http://docs.juliadiffeq.org/stable/tutorials/sde_example.html#Example-4:-Systems-of-SDEs-with-Non-Diagonal-Noise-1
%using my stochastic integrator
    npaths=10000;
    nnoise=4;
    u0=ones(2,1);
    ti=0.0;
    tf=1.0;
    dt=(tf-ti)/10000;
    
    S=zeros(npaths,2,101);
    T=zeros(npaths,101);
    
    parfor m=1:npaths
        [S1,T1]=rk4int(u0,@f,@g,nnoise,ti,tf,dt,{@sample},[100]);
        S(m,:,:)=S1{1};
        T(m,:)=T1{1};
    end
    
    

end

function G=g(u,t)
    G=[0.3*u(1) 0.6*u(1) 0.9*u(1) 0.12*u(2);
        1.2*u(1) 0.2*u(2) 0.3*u(2) 1.8*u(2)];
end

function F=f(u,t)
    F=1.01*u;
end

function s=sample(z,t)
    s=z;
end