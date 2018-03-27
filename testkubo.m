%test simple kubo oscillator with rk4 integrator

function [S,T]= testkubo()
    z0=ones(10000,1);
    ti=0;
    tf=10;
    dt=(tf-ti)/1000;
    
    function F=f(z,t)
        F=0;
    end

    function G=g(z,t)
        G=diag(1j*z);
    end
    
    
    function s=sample(z,t)
        s=z;
    end

    
    [S,T]=rk4int(z0,@f,@g,ti,tf,dt,{@sample},[100]);
    
    
end