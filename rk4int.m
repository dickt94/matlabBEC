%rk4int - stochastic RK4 integrator.
%integrates the Stratonovich equation dc=f(c,t)dt+g(c,t)dW from ti to tf,
%using timestep dt and an RK4 algorithm.

%arguments:
%ci - initial field, vector of size N
%f - function handle, takes arguments c,t, returns vector of size N
%nproc - number of noise processes
%g - nproc length cell array of function handles, takes arguments c,t, returns N*nnoise matrix or
%nnoise*1 vector if diag is 1.
%nnoise - sizes of noise processes - vector of length nproc
%seed - nproc length vector of seeds for each noise process. Defaults to
%'shuffle' seed if set to zero.
%ti - initial time
%tf - final time
%nsteps - number of steps to take between ti and tf
%moments - cell array function handles to be sampled, size M. Functions
%accept c,t as arguments.
%samples - array of numbers of samples to be taken
%diag - array of bool, is each of the noises diagonal?
%filters - cell array of functions that accept c,t as arguments
%and return a new c. These are executed every step at step end, in 
%sequence. The field is replaced by that returned by the filter.
%ipop - interaction picture operator, N*N matrix if ipdiag is 0, else N*1
%vector
%ipdiag - is the ip operator diagonal?


%also performs a half-step integration with timestep dt/2 to check
%convergence.


%outputs:
%sdata - cell array containing sample results. The nth element of sdata
%corresponds to the nth function passed in moments, and should be an array
%of dimensionality one greater than that of the array returned by the
%sample function. These arrays are formed by concatenating the sample
%results along their first singleton dimension.
%tdata - cell array containing 1D arrays of sampling times. nth element
%corresponds to nth element of moments.

%Richard Taylor, 2018.

%todo: move sampling to half-step to give more precision I guess.

function [sdata,tdata]=rk4int(ci,f,ipop,ipdiag,nproc,g,nnoise,seed,ti,tf,nsteps,moments,samples,diag,filters)
    
    %create random streams with seeds specified. Calling from any one
    %stream should not affect the other.
    streams=cell(nproc,1);
    for count=1:nproc
        %create a random stream with the specified seed
        streams{count}=RandStream.create('dsfmt19937', 'NumStreams',1,'Seed',seed(count));
    end
    
    function v=G(d,u,dX)
        v=f(d,u);
        for cnt=1:nproc
            if diag(cnt)
                v=v+g{cnt}(d,u).*dX{cnt};
            else
                v=v+g{cnt}(d,u)*dX{cnt};
            end
        end
        
    end
    
    c=ci;
    t=ti;
    
    %figure out when to take samples
    stepcount=0;
    dt=(tf-ti)/nsteps;
    sint=round(nsteps./samples);
    
    if ~isempty(ipop)
        if ipdiag
            evolip=exp(ipop*dt/2);
        else
            evolip=expm(ipop*dt/2);
        end
    else
        evolip=1;
    end
    
    %initialise sample return
    sdataf=cell(length(moments),1);
    tdataf=cell(length(moments),1);
    sampdim=cell(length(moments),1);
    
    %initial sample
    for i=1:length(moments)
        %sampling concatenates along the first singleton dimension, or if
        %there is none, makes an extra dimension to concatenate along. This
        %should result in correct handling of sample functions that return
        %column vectors and row vectors, or higher-dimensional arrays.
        samp=moments{i}(c,t);
        ssamp=size(samp);
        sampdim{i}=1;
        for p=1:length(ssamp)
            if ssamp(p)==1
                break;
            end
            sampdim{i}=sampdim{i}+1;
        end
        sdataf{i}=cat(sampdim{i},sdataf{i},samp);
        tdataf{i}(length(tdataf{i})+1)=t;
        fprintf("Sampled moment %d at t=%f\n",i,t);
    end
    
    sigma=sqrt(2/dt);
    

    while stepcount < nsteps
        dW=cell(nproc,1);
        for count=1:nproc
            dW1=sigma*randn(streams{count},[nnoise(count) 1]);
            dW2=sigma*randn(streams{count},[nnoise(count) 1]);%roll twice - this keeps the seeding the same for halfstep
        
            dW{count}=(dW1+dW2)/2;
        end


    
        if ipdiag
            cI=evolip.*c;
            ck=evolip.*G(c,t,dW)*dt;%k1
            c=cI+ck/6;
            ck=ck/2+cI;
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k2
            c=c+ck/3;
            ck=ck/2+cI;
            ck=G(ck,t,dW)*dt;%k3
            c=c+ck/3;
            ck=evolip.*(ck+cI);
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k4
            c=evolip.*c+ck/6;
        else
            cI=evolip*c;
            ck=evolip*G(c,t,dW)*dt;%k1
            c=cI+ck/6;
            ck=ck/2+cI;
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k2
            c=c+ck/3;
            ck=ck/2+cI;
            ck=G(ck,t,dW)*dt;%k3
            c=c+ck/3;
            ck=evolip*(ck+cI);
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k4
            c=evolip*c+ck/6;
        end
        
        %do filters
        for i=1:length(filters)
            c=filters{i}(c,t);
        end
        
        
        stepcount=stepcount+1;
        takeSample=~mod(stepcount,sint);
        for i=1:length(takeSample)
            if takeSample(i)
                sdataf{i}=cat(sampdim{i},sdataf{i},moments{i}(c,t));
                tdataf{i}(length(tdataf{i})+1)=t;
                fprintf("Sampled moment %d at t=%f\n",i,t);
            end
        end

    end

    
    %half-step

    fprintf("Beginning half-step integration...\n");
    
    c=ci;
    t=ti;
    dt=dt/2;
    
    %recalc the IP matrix for halfstep
    if ~isempty(ipop)
        if ipdiag
            evolip=exp(ipop*dt/2);
        else
            evolip=expm(ipop*dt/2);
        end
    else
        evolip=1;
    end
    
    %same noises as for full-step
    for count=1:nproc
        %create a random stream with the specified seed
        streams{count}=RandStream.create('dsfmt19937', 'NumStreams',1,'Seed',seed(count));
    end
    
    
    %initialise half-step sample return
    sdata=cell(length(moments),1);
    tdata=cell(length(moments),1);
    sampdim=cell(length(moments),1);
    
    %initial sample
    for i=1:length(moments)
        %sampling concatenates along the first singleton dimension, or if
        %there is none, makes an extra dimension to concatenate along. This
        %should result in correct handling of sample functions that return
        %column vectors and row vectors, or higher-dimensional arrays.
        samp=moments{i}(c,t);
        ssamp=size(samp);
        sampdim{i}=1;
        for p=1:length(ssamp)
            if ssamp(p)==1
                break;
            end
            sampdim{i}=sampdim{i}+1;
        end
        sdata{i}=cat(sampdim{i},sdata{i},samp);
        tdata{i}(length(tdata{i})+1)=t;
        fprintf("Sampled moment %d at t=%f\n",i,t);
    end
    
    serr=zeros([length(moments) 1]);
    
    
    stepcount=0;
    while stepcount<nsteps*2
            
        dW=cell(nproc,1);
        for count=1:nproc
            dW{count}=sigma*randn(streams{count},[nnoise(count) 1]);
        end
    

    
        if ipdiag
            cI=evolip.*c;
            ck=evolip.*G(c,t,dW)*dt;%k1
            c=cI+ck/6;
            ck=ck/2+cI;
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k2
            c=c+ck/3;
            ck=ck/2+cI;
            ck=G(ck,t,dW)*dt;%k3
            c=c+ck/3;
            ck=evolip.*(ck+cI);
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k4
            c=evolip.*c+ck/6;
        else
            cI=evolip*c;
            ck=evolip*G(c,t,dW)*dt;%k1
            c=cI+ck/6;
            ck=ck/2+cI;
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k2
            c=c+ck/3;
            ck=ck/2+cI;
            ck=G(ck,t,dW)*dt;%k3
            c=c+ck/3;
            ck=evolip*(ck+cI);
            t=t+dt/2;
            ck=G(ck,t,dW)*dt;%k4
            c=evolip*c+ck/6;
        end
        
        
        %do filters
        for i=1:length(filters)
            c=filters{i}(c,t);
        end
        
        
        stepcount=stepcount+1;
        takeSample=~mod(stepcount,sint*2);%half-step samples
        for i=1:length(takeSample)
            if takeSample(i)
                samp=moments{i}(c,t);
                
                sdata{i}=cat(sampdim{i},sdata{i},samp);
                tdata{i}(length(tdata{i})+1)=t;
                
                %todo calculate error for ith moment and check if bigger
                %than maximum
                
                szfsamp=size(sdataf{i});
                inds=cell(1,length(szfsamp));
                for p=1:length(szfsamp)
                    if p==sampdim{i}
                        inds{p}=round(stepcount/(sint(i)*2))+1;
                    else
                        inds{p}=1:szfsamp(p);
                    end
                end
                fullsamp=sdataf{i}(inds{:});
                steperr=norm(fullsamp-samp,inf);
                if steperr>serr(i)
                    serr(i)=steperr;
                end
                fprintf("Sampled moment %d at t=%f\n",i,t);
            end
        end
        
    end
    
    for p=1:length(moments)
        fprintf("Maximum step error in moment %d was %f\n",p,serr(p));
    end
    
end

