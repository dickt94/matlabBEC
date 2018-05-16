%test an NPW representation of the harmonic oscillator under continuous
%number measurement. Should give the same results as the direct integration
%of this problem.

%evolves the normalised weight equations.

function [S,T]=testsingle()
    %set up constants
    nbar=10;
    phibar=0;
    gamma=1.0;
    sgam=sqrt(gamma);

    %k samples of n - drawn from Poissonian distribution to represent
    %sampling an initial coherent state
    k=10000;
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
        w=exp(field(2*k+1:3*k));
        expectn=(sum(w.*n)/sum(w));
        
        %phi=field(k+1:2k);
        %logw=field(2*k+1:3*k);
        F=[zeros(k,1);zeros(k,1);2*sgam*(n-expectn)*realnoise(t)+gamma*((-2*n+4*expectn).*n-(-2*(sum(w.*n.*n)/sum(w))+4*expectn^2))];
    end

    function G=g(field,t)
        G=[zeros(k,1);sgam*ones(k,1);zeros(k,1)];
    end

    function nb=meann(field,t)
        w=exp(field(2*k+1:3*k));
        nb=sum(field(1:k).*w)/sum(w);
    end

    function vn=varn(field,t)
        w=exp(field(2*k+1:3*k));
        n=field(1:k);
        %bootstrap
        nsamp=100;
        slen=floor(k/nsamp);
        means=zeros([nsamp 1]);
        for samp=1:nsamp
            means(samp)=sum(n((samp-1)*slen+1:samp*slen).*w((samp-1)*slen+1:samp*slen))/sum(w((samp-1)*slen+1:samp*slen));
        end
        vn=std(means)/sqrt(nsamp-1);
    end

    function fs=fieldsamp(field,t)
        fs=field;
    end

    function ws=weightsamp(field,t)
        ws=[max(field(2*k+1:3*k)),min(field(2*k+1:3*k))];
    end

    %maximum log difference (ln(wmax/wmin)) between smallest and largest
    %weights
    breedgap=100;
    breedcount=0;
    function field=breed(field,t)
        [minw,mini]=min(field(2*k+1:3*k));
        [maxw,maxi]=max(field(2*k+1:3*k));
        
        %normalise
        %field(2*k+1:3*k)=field(2*k+1:3*k)-maxw;
        
        %breed
        while maxw-minw>breedgap
            field(2*k+mini)=maxw-log(2);
            field(2*k+maxi)=maxw-log(2);
            
            %copy n and phi across
            field(mini)=field(maxi);
            field(k+mini)=field(k+maxi);
            
            breedcount=breedcount+1
            
            [minw,mini]=min(field(2*k+1:3*k));
            [maxw,maxi]=max(field(2*k+1:3*k));
        end
        
        %no normalisation required because the weight equations are already
        %normalised
    end


    %generate a timestep-normalised real noise for the fictitious path
    %evolution
    interval=6;
    nsteps=15000;
    dt=interval/nsteps;
    sig=sqrt(2/dt);
    
    %make all the real noises now so as not to be interfered with by the
    %rk4int 
    rng(314159265);
    rnoises=normrnd(0,sig,[2*(nsteps+1) 1]);
    rnoisesf=zeros([nsteps+1 1]);
    for count=0:nsteps
        rnoisesf(count+1)=(rnoises(2*count+1)+rnoises(2*count+2))/2;
    end
    
    function dW=realnoise(t)

        stepcount=round(t/dt)+1;
        dW=rnoisesf(stepcount);

    end


    [S,T]=rk4int(initf,@f,@g,nnoise,0,interval,nsteps,{@meann,@varn,@fieldsamp,@weightsamp},[100 100 100 100],true,{@breed});

end
        
        