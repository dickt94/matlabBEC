%test simple kubo oscillator with rk4 integrator

function [S,T]= testkubo()
    npaths=100;

    S=cell(npaths,1);
    T=cell(npaths,1);
    
    z0=1;
    ti=0;
    tf=10;
    dt=(tf-ti)/1000;
    


    parfor m=1:npaths
        [S{m},T{m}]=rk4int(z0,@f,@g,1,ti,tf,dt,{@sample},[100]);
    end
    
    
end

function F=f(z,t)
    F=0;
end

function G=g(z,t)
    G=1j*z;
end
    
    
function s=sample(z,t)
    s=z;
end