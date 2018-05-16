%rk4int_direct - stochastic RK4 integrator. Special version for comparison
%with NPW simulation.
%integrates the Stratonovich equation dc=f(c,t)dt+g(c,t)dW from ti to tf,
%using timestep dt and an RK4 algorithm.

%arguments:
%ci - initial field, vector of size N
%f - function handle, takes arguments c,t, returns vector of size N
%g - function handle, takes arguments c,t, returns N*nnoise matrix or
%nnoise*1 vector if diag is 1.
%nnoise - number of noise processes
%ti - initial time
%tf - final time
%nsteps - number of steps to take between ti and tf
%moments - cell array function handles to be sampled, size M. Functions
%accept c,t as arguments.
%samples - array of numbers of samples to be taken
%diag - bool, are the noises diagonal?
%filters - cell array of functions that accept c,t as arguments and return 
%a new c. These are executed every step at step end, in sequence. The field
%is replaced by that returned by the filter.

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

function [sdata,tdata]=rk4int(ci,f,g,nnoise,ti,tf,nsteps,moments,samples,diag,filters)
    
    function v=G(d,u,dX)
        if diag
            v=(f(d,u)+g(d,u).*dX);
        else
            v=(f(d,u)+g(d,u)*dX);
        end
        
    end
    
    c=ci;
    t=ti;
    
    
    
    %figure out when to take samples
    stepcount=0;
    dt=(tf-ti)/nsteps;
    sint=round(nsteps./samples);
    
    %initialise sample return
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
        fprintf("Sampled function %d at t=%f\n",i,t);
    end
    
    sigma=sqrt(2/dt);
    
    %roll noises now with fixed seed, for compatibility with NPW
    %simulation.
    rng(314159265);
    rnoises=normrnd(0,sigma,[2*(nsteps+1) 1]);
    

    while stepcount < nsteps
        dW1=rnoises(2*stepcount+1);
        dW2=rnoises(2*stepcount+2);%roll twice - this keeps the seeding the same for halfstep
        
        dW=(dW1+dW2)/2;


    
        cI=c;
        ck=G(c,t,dW)*dt;%k1
        c=c+ck/6;
        ck=ck/2+cI;
        t=t+dt/2;
        ck=G(ck,t,dW)*dt;%k2
        c=c+ck/3;
        ck=ck/2+cI;
        ck=G(ck,t,dW)*dt;%k3
        c=c+ck/3;
        ck=ck+cI;
        t=t+dt/2;
        ck=G(ck,t,dW)*dt;%k4
        c=c+ck/6;
        
        %do filters
        for i=1:length(filters)
            c=filters{i}(c,t);
        end
        
        
        stepcount=stepcount+1;
        takeSample=~mod(stepcount,sint);
        for i=1:length(takeSample)
            if takeSample(i)
                sdata{i}=cat(sampdim{i},sdata{i},moments{i}(c,t));
                tdata{i}(length(tdata{i})+1)=t;
                fprintf("Sampled function %d at t=%f\n",i,t);
            end
        end

    end
    cfull=c;

    
    %half-step

    c=ci;
    t=ti;
    dt=dt/2;
    
    
    stepcount=0;
    while stepcount<nsteps*2
            
        dW=rnoises(stepcount+1);
    

    
        cI=c;
        ck=G(c,t,dW)*dt;%k1
        c=c+ck/6;
        ck=ck/2+cI;
        t=t+dt/2;
        ck=G(ck,t,dW)*dt;%k2
        c=c+ck/3;
        ck=ck/2+cI;
        ck=G(ck,t,dW)*dt;%k3
        c=c+ck/3;
        ck=ck+cI;
        t=t+dt/2;
        ck=G(ck,t,dW)*dt;%k4
        c=c+ck/6;
        
        
        %do filters
        for i=1:length(filters)
            c=filters{i}(c,t);
        end
        
        
        stepcount=stepcount+1;
        
    end
    
    norm(c-cfull)
    
end

