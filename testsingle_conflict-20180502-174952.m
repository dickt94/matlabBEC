%test an NPW representation of the harmonic oscillator under continuous
%number measurement. Should give the same results as the direct integration
%of this problem.

function [S,T]=testsingle()
    %set up constants
    nbar=10;
    phibar=0;
    gamma=1.0;
    sgam=sqrt(gamma);

    
    %k samples of n - drawn from Poissonian distribution to represent
    %sampling an initial coherent state
    k=100;
    nnoise=3*k;%k fictitious noises but it's easier to just roll them for the whole vector
    rng('shuffle');
    for c=1:k
        initn(c)=poissrnd(nbar);
        if initn(c)==0
            initp(c)=2*pi*rand;
        else
            initp(c)=normrnd(phibar,sqrt(psi(1,c+1)));
        end
    end
    initn=initn.';
    initp=initp.';
    
    initW=ones(k,1);
    
    %weights integrated in log space
    initw=log(initW);
    
    initf=[initn;initp;initw];
    
    function F=f(field,t)
        n=field(1:k);
        %phi=field(k+1:2k);
        %logw=field(2*k+1:3*k);
        F=[zeros(k,1);zeros(k,1);gamma*(-2*n.*n+4*mean(n.*exp(field(2*k+1:3*k)))*n)+2*sgam*n*realnoise(t)];
    end

    function G=g(field,t)
        G=[zeros(k,1);sgam*ones(k,1);zeros(k,1)];
    end

    function nb=meann(field,t)
        nb=mean(field(1:k).*field(2*k+1:3*k));
    end


    %generate a timestep-normalised real noise for the fictitious path
    %evolution
    interval=6;
    nsteps=20000;
    dt=interval/nsteps;
    sig=sqrt(2/dt);
    
    %make all the real noises now so as not to be interfered with by the
    %rk4int 
    rng(314159265);
    rnoises=normrnd(0,sig,[2*nsteps 1]);
    rnoisesf=zeros([nsteps 1]);
    for count=0:nsteps-1
        rnoisesf(count+1)=(rnoises(2*count+1)+rnoises(2*count+2))/2;
    end
    
    function dW=realnoise(t)
            stepcount=round(t/dt)
            dW=rnoisesf(stepcount);
    end


    [S,T]=rk4int(initf,@f,@g,nnoise,0,interval,nsteps,{@meann},[100],true);

end
        
        