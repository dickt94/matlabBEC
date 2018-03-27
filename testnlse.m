%test the RK4 integrator with a nonlinear schroedinger equation
%isolates the deterministic part of the evolution

function [samples,times]=testnlse()
    energy=4;
    vel=0.3;
    hwhm=1.0;
    
    tau=linspace(-6,6,129);
    tau=tau(1:end-1);
    ktau=ifftshift((-64:63)*2*pi/12)%k-space grid (hopefully correct)
    
    
    %initial wavefunction
    w0=hwhm*sqrt(2/log(2));
    amp=sqrt(energy/w0/sqrt(pi/2));
    phi0=amp*exp(-tau.*tau/w0/w0).*exp(1i*vel*tau);
    
    %damping vector
    
    Gamma=1.0*(1-exp(-(tau.*tau/4.0/4.0).^10));
    
    function F=f(phi,t)
        phik=fft(phi);
        deriv2term=ifft(-0.5i*ktau.^2.*phik);
        
        F=deriv2term-phi.*Gamma+1i.*conj(phi).*phi.*phi;
    end

    function G=g(phi,t)
        G=0;
    end

    function d=dens(phi,t)
        d=conj(phi).*phi;
    end

    function n=norm(phi,t)
        dx=12/128;
        n=sum(conj(phi).*phi)*dx;
    end


    function d2t=kterm(phi,t)
        phik=fft(phi);
        d2t=ifft(-0.5i*ktau.^2.*phik);
    end


    ti=0;
    tf=20;
    dt=(tf-ti)/10000;
    [samples,times]=rk4int(phi0,@f,@g,ti,tf,dt,{@dens,@norm,@kterm},[100 100 100]);


end