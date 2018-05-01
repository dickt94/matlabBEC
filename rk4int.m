%rk4int - stochastic RK4 integrator.
%integrates the Stratonovich equation dc=f(c,t)dt+g(c,t)dW from ti to tf,
%using timestep dt and an RK4 algorithm.

%arguments:
%ci - initial field, vector of size N
%f - function handle, takes arguments c,t, returns vector of size N
%g - function handle, takes arguments c,t, returns N*nnoise matrix
%nnoise - number of noise processes
%ti - initial time
%tf - final time
%dt - desired time step
%moments - cell array function handles to be sampled, size M. Functions
%accept c,t as arguments.
%samples - array of numbers of samples to be taken

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

function [sdata,tdata]=rk4int(ci,f,g,nnoise,ti,tf,dt,moments,samples)
    
    function v=G(d,u,dX)
        v=(f(d,u)+g(d,u)*dX);
    end
    
    c=ci;
    t=ti;
    
    
    rng('shuffle');
    seed=rng;
    
    %figure out when to take samples
    stepcount=0;
    nsteps=(tf-ti)/dt;
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
    
    while 1
        if dt < tf-t
            
            dW1=normrnd(0,sigma,[nnoise 1]);
            dW2=normrnd(0,sigma,[nnoise 1]);%roll twice - this keeps the seeding the same for halfstep
            
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
            
        else
            dt1=tf-t;
            dW1=normrnd(0,1,[nnoise 1]);
            dW2=normrnd(0,1,[nnoise 1]);%roll twice - this keeps the seeding the same for halfstep
            
            dW=(dW1+dW2)/sqrt(2*dt1);
    

    
            cI=c;
            ck=G(c,t,dW)*dt1;%k1
            c=c+ck/6;
            ck=ck/2+cI;
            t=t+dt1/2;
            ck=G(ck,t,dW)*dt1;%k2
            c=c+ck/3;
            ck=ck/2+cI;
            ck=G(ck,t,dW)*dt1;%k3
            c=c+ck/3;
            ck=ck+cI;
            t=t+dt1/2;
            ck=G(ck,t,dW)*dt1;%k4
            c=c+ck/6;
            break;
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
    rng(seed);%same noises as for full-step
    while 1
        if dt < tf-t
            
            dW=normrnd(0,sigma,[nnoise 1]);
    

    
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
            
        else
            dt1=tf-t;
            dW=normrnd(0,1,[nnoise 1])/sqrt(dt1);
    

    
            cI=c;
            ck=G(c,t,dW)*dt1;%k1
            c=c+ck/6;
            ck=ck/2+cI;
            t=t+dt1/2;
            ck=G(ck,t,dW)*dt1;%k2
            c=c+ck/3;
            ck=ck/2+cI;
            ck=G(ck,t,dW)*dt1;%k3
            c=c+ck/3;
            ck=ck+cI;
            t=t+dt1/2;
            ck=G(ck,t,dW)*dt1;%k4
            c=c+ck/6;
            break;
        end
    end
    
    norm(c-cfull)
    
end

