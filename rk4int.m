%rk4int - stochastic RK4 integrator.
%integrates the Stratonovich equation dc=f(c,t)dt+g(c,t)dW from ti to tf,
%using timestep dt and an RK4 algorithm.

%arguments:
%ci - initial field, column vector of size N
%f - function handle, takes arguments c,t, returns column vector of size N
%g - function handle, takes arguments c,t, returns N*N matrix
%ti - initial time
%tf - final time
%dt - desired time step
%moments - cell array function handles to be sampled, size M. Functions
%accept c,t as arguments.
%samples - array of numbers of samples to be taken

%also performs a half-step integration with timestep dt/2 to check
%convergence.

%Richard Taylor, 2018.

function [sdata,tdata]=rk4int(ci,f,g,ti,tf,dt,moments,samples)
    
    function v=G(d,u,dX)
        v=(f(d,u)+g(d,u)*dX);
    end
    
    c=ci;
    t=ti;
    
    
    rng(0);
    
    %figure out when to take samples
    stepcount=0;
    nsteps=(tf-ti)/dt;
    sint=round(nsteps./samples);
    
    %initialise sample return
    sdata=cell(length(moments),1);
    tdata=cell(length(moments),1);
    
    %initial sample
    for i=1:length(moments)
        sdata{i}=cat(1,sdata{i},moments{i}(c,t));
        tdata{i}(length(tdata{i})+1)=t;
    end
    
    while 1
        if dt < tf-t
            
            dW1=normrnd(0,1,size(c));
            dW2=normrnd(0,1,size(c));%roll twice - this keeps the seeding the same for halfstep
            
            dW=(dW1+dW2)/sqrt(2*dt);
    

    
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
            dW1=normrnd(0,1,size(c));
            dW2=normrnd(0,1,size(c));%roll twice - this keeps the seeding the same for halfstep
            
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
                sdata{i}=cat(2,sdata{i},moments{i}(c,t));
                tdata{i}(length(tdata{i})+1)=t;
            end
        end

    end
    cfull=c;

    
    %half-step
    c=ci;
    t=ti;
    dt=dt/2;
    rng(0);%same noises as for full-step
    while 1
        if dt < tf-t
            
            dW=normrnd(0,1,size(c))/sqrt(dt);
    

    
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
            dW=normrnd(0,1,size(c))/sqrt(dt1);
    

    
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

