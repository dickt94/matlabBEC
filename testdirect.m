%single-mode harmonic oscillator under number measurement, useful for
%comparison to NPW.

function [S,T]=testdirect()
    nmodes=30;
    nnoise=1;
    nbar=10;
    
    alph=sqrt(nbar);
    gamma=1.0;
    sgam=sqrt(gamma);
    
    %define field operators
    a=diag(sqrt(1:nmodes-1),1);
    ad=conj(a.');
    ada=ad*a;
    
    %coherent state initial
    n=0:nmodes-1;
    psi0=exp(-abs(alph)^2/2)*(alph.^n)./sqrt(factorial(n));
    
    rho0=conj(psi0.')*psi0;
    
    %vector for integration
    rhovec0=rho0(:);
    
    function F=f(rhovec,t)
        rho=reshape(rhovec,[nmodes nmodes]);
        
        expectn=trace(rho*ada);
        expectn2=trace(rho*ada*ada);
        
        Hrho=ada*rho+rho*ada-2*expectn*rho;
        H2rho=ada*ada*rho+rho*ada*ada-2*expectn2*rho;
        
        %gamma*([D[ad*a]rho+C[ad*a]rho)
        Fmat=gamma*(ada*rho*ada-1/2 * (ada*ada*rho+rho*ada*ada) + 2*expectn*Hrho-1/2*H2rho+expectn2*rho-ada*rho*ada);
        
        F=Fmat(:);
        
        
    end

    function G=g(rhovec,t)
        rho=reshape(rhovec,[nmodes nmodes]);
        %sqrt(gamma)*H[ad*a]rho
        Gmat=sgam*(ada*rho+rho*ada-trace(rho*2*ada)*rho);
        G=Gmat(:);
    end

    interval=6;
    nsteps=15000;
    
    function sample=s(rhovec,t)
        rho=reshape(rhovec,[nmodes nmodes]);
        sample=trace(rho*ada);
    end
    function sample=s2(rhovec,t)
        sample=rhovec;
    end

    function sample=stratcorrection(rhovec,t)
        rho=reshape(rhovec,[nmodes nmodes]);
        expectn=trace(rho*ada);
        expectn2=trace(rho*ada*ada);
        
        Hrho=ada*rho+rho*ada-2*expectn*rho;
        H2rho=ada*ada*rho+rho*ada*ada-2*expectn2*rho;
        
        samplemat=gamma*(2*expectn*Hrho-1/2*H2rho+expectn2*rho-ada*rho*ada);
        sample=samplemat(:);
    end

    function sample=innov(rhovec,t)
        rho=reshape(rhovec,[nmodes nmodes]);
        samplemat=sgam*(ada*rho+rho*ada-trace(rho*2*ada)*rho);
        sample=samplemat(:);
    end
    
    [S,T]=rk4int_direct(rhovec0,@f,@g,nnoise,0,interval,nsteps,{@s,@s2,@stratcorrection,@innov},[100 100 100 100],false,{});
    
end