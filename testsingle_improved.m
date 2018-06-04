%test an NPW representation of the harmonic oscillator under continuous
%number measurement. Should give the same results as the direct integration
%of this problem.

function [S,T]=testsingle_improved()
    %set up constants
    nbar=10;
    phibar=0;
    gamma=1.0;
    sgam=sqrt(gamma);

    %k samples of n - drawn from Poissonian distribution to represent
    %sampling an initial coherent state
    k=15000;
    nnoise=[3*k 1];%k fictitious noises but it's easier to just roll them for the whole vector
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
        F=[zeros(k,1);zeros(k,1);gamma*(-2*n+4*expectn).*n];
    end

    function G=g1(field,t)
        G=[zeros(k,1);sgam*ones(k,1);zeros(k,1)];
    end

    function G=g2(field,t)
        G=[zeros(k,1);zeros(k,1);2*sgam*field(1:k)];
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
        
        %normalise
        maxw=max(field(2*k+1:3*k));
        field(2*k+1:3*k)=field(2*k+1:3*k)-maxw;
        
    end


    %integration interval and number of steps to take
    interval=6;
    nsteps=15000;



    [S,T]=rk4int(initf,@f,2,{@g1 @g2},nnoise,[231472190 314159265],0,interval,nsteps,{@meann,@varn,@fieldsamp,@weightsamp},[100 100 100 100],[true true],{@breed});

end
        
        