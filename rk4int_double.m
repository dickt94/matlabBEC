%rk4int - stochastic RK4 integrator.
%integrates the Stratonovich equation dc=f(c,t)dt+g(c,t)dW from ti to tf,
%using timestep dt and an RK4 algorithm.

%arguments:
%ci - initial field, size N
%f - function handle, takes arguments c,t, returns vector of size N
%nproc - number of noise processes
%g - nproc length cell array of function handles, takes arguments c,t, dW
%and returns array of size N
%nnoise - sizes of noise processes - cell array of matrices with noise
%sizes in each dimension
%seed - nproc length vector of seeds for each noise process. Defaults to
%'shuffle' seed if set to zero.
%ti - initial time
%tf - final time
%nsteps - number of steps to take between ti and tf
%moments - cell array function handles to be sampled, size M. Functions
%accept c,t as arguments.
%samples - array of numbers of samples to be taken
%filters - cell array of functions that accept c,t as arguments
%and return a new c. These are executed every step at step end, in 
%sequence. The field is replaced by that returned by the filter.
%ipop - interaction picture operator, N*N matrix if ipdiag is 0, else N*1
%vector
%ipdiag - is the ip operator diagonal?
%errorcheck - bool - perform half-step integration?

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
%this is a hacky double-vector version for system-filter separation.

%todo: move sampling to half-step to give more precision I guess.

function [sdata,tdata]=rk4int_double(ci1,ci2,f1,f2,ipop1,ipop2,ipdiag1,ipdiag2,nproc,g1,g2,nnoise,seed,ti,tf,nsteps,moments,samples,filters,errorcheck)
    
    %create random streams with seeds specified. Calling from any one
    %stream should not affect the other.
    streams=cell(nproc,1);
    for count=1:nproc
        %create a random stream with the specified seed
        streams{count}=RandStream.create('dsfmt19937', 'NumStreams',1,'Seed',seed(count));
    end
    
    function v=G1(d1,d2,u,dX)
        v=f1(d1,d2,u);
        for cnt=1:nproc
            v=v+g1{cnt}(d1,d2,u,dX{cnt});
        end
        
    end
    
    function v=G2(d1,d2,u,dX)
        v=f2(d1,d2,u);
        for cnt=1:nproc
            v=v+g2{cnt}(d1,d2,u,dX{cnt});
        end
        
    end

    c1=ci1;
    c2=ci2;
    t=ti;
    
    %figure out when to take samples
    stepcount=0;
    dt=(tf-ti)/nsteps;
    sint=round(nsteps./samples);
    
    if ~isempty(ipop1)
        if ipdiag1
            evolip1=exp(ipop1*dt/2);
        else
            evolip1=expm(ipop1*dt/2);
        end
    else
        evolip1=1;
    end
    if ~isempty(ipop2)
        if ipdiag2
            evolip2=exp(ipop1*dt/2);
        else
            evolip2=expm(ipop1*dt/2);
        end
    else
        evolip2=1;
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
        samp=moments{i}(c1,c2,t);
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
    
    %initialise the noise sampling
    %sdataf{length(moments)+1}=[];
    %tdataf{length(moments)+1}=[];
    
    sigma=sqrt(2/dt);
    

    while stepcount < nsteps
        dW=cell(nproc,1);
        for count=1:nproc
            dW1=sigma*randn(streams{count},[nnoise{count} 1]);
            dW2=sigma*randn(streams{count},[nnoise{count} 1]);%roll twice - this keeps the seeding the same for halfstep
        
            dW{count}=(dW1+dW2)/2;
        end


        %system-filter separation - integrate system and filter (c1 and c2)
        %in parallel
    
        if ipdiag1 && ipdiag2
            cI1=evolip1.*c1;
            cI2=evolip2.*c2;
            ck1=evolip1.*G1(c1,c2,t,dW)*dt;%k1
            ck2=evolip2.*G2(c1,c2,t,dW)*dt;%k1
            c1=cI1+ck1/6;
            c2=cI2+ck2/6;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k2
            ck2=G2(ck1,ck2,t,dW)*dt;%k2
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k3
            ck2=G2(ck1,ck2,t,dW)*dt;%k3
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=evolip1.*(ck1+cI1);
            ck2=evolip2.*(ck2+cI2);
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k4
            ck2=G2(ck1,ck2,t,dW)*dt;%k4
            c1=evolip1.*c1+ck1/6;
            c2=evolip2.*c2+ck2/6;
        elseif ipdiag1
            cI1=evolip1.*c1;
            cI2=evolip2*c2;
            ck1=evolip1.*G1(c1,c2,t,dW)*dt;%k1
            ck2=evolip2*G2(c1,c2,t,dW)*dt;%k1
            c1=cI1+ck1/6;
            c2=cI2+ck2/6;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k2
            ck2=G2(ck1,ck2,t,dW)*dt;%k2
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k3
            ck2=G2(ck1,ck2,t,dW)*dt;%k3
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=evolip1.*(ck1+cI1);
            ck2=evolip2*(ck2+cI2);
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k4
            ck2=G2(ck1,ck2,t,dW)*dt;%k4
            c1=evolip1.*c1+ck1/6;
            c2=evolip2*c2+ck2/6;
        elseif ipdiag2
            cI1=evolip1*c1;
            cI2=evolip2.*c2;
            ck1=evolip1*G1(c1,c2,t,dW)*dt;%k1
            ck2=evolip2.*G2(c1,c2,t,dW)*dt;%k1
            c1=cI1+ck1/6;
            c2=cI2+ck2/6;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k2
            ck2=G2(ck1,ck2,t,dW)*dt;%k2
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k3
            ck2=G2(ck1,ck2,t,dW)*dt;%k3
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=evolip1*(ck1+cI1);
            ck2=evolip2.*(ck2+cI2);
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k4
            ck2=G2(ck1,ck2,t,dW)*dt;%k4
            c1=evolip1*c1+ck1/6;
            c2=evolip2.*c2+ck2/6;            
        else
            cI1=evolip1*c1;
            cI2=evolip2*c2;
            ck1=evolip1*G1(c1,c2,t,dW)*dt;%k1
            ck2=evolip2*G2(c1,c2,t,dW)*dt;%k1
            c1=cI1+ck1/6;
            c2=cI2+ck2/6;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k2
            ck2=G2(ck1,ck2,t,dW)*dt;%k2
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=ck1/2+cI1;
            ck2=ck2/2+cI2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k3
            ck2=G2(ck1,ck2,t,dW)*dt;%k3
            c1=c1+ck1/3;
            c2=c2+ck2/3;
            ck1=evolip1*(ck1+cI1);
            ck2=evolip2*(ck2+cI2);
            t=t+dt/2;
            ck1=G1(ck1,ck2,t,dW)*dt;%k4
            ck2=G2(ck1,ck2,t,dW)*dt;%k4
            c1=evolip1*c1+ck1/6;
            c2=evolip2*c2+ck2/6;
        end
        
        %do filters
        for i=1:length(filters)
            [c1,c2]=filters{i}(c1,c2,t);
        end
        
        
        stepcount=stepcount+1;
        takeSample=~mod(stepcount,sint);
        for i=1:length(takeSample)
            if takeSample(i)
                sdataf{i}=cat(sampdim{i},sdataf{i},moments{i}(c1,c2,t));
                tdataf{i}(length(tdataf{i})+1)=t;
                fprintf("Sampled moment %d at t=%f\n",i,t);
            end
        end

                
        %sample the noise processes
        %takenSample=~mod(stepcount,nsteps/noisesamps);
        %if takenSample
        %    sdataf{length(moments)+1}=[sdataf{length(moments)+1} dW];
        %    tdataf{length(moments)+1}=[tdataf{length(moments)+1} t];
        %    fprintf("sampled noise\n");
        %end
        
    end

    
    %half-step
    if errorcheck
        fprintf("Beginning half-step integration...\n");

        c1=ci1;
        c2=ci2;
        t=ti;
        dt=dt/2;

        %recalc the IP matrix for halfstep
        if ~isempty(ipop1)
            if ipdiag1
                evolip1=exp(ipop1*dt/2);
            else
                evolip1=expm(ipop1*dt/2);
            end
        else
            evolip1=1;
        end
        if ~isempty(ipop2)
            if ipdiag2
                evolip2=exp(ipop2*dt/2);
            else
                evolip2=expm(ipop2*dt/2);
            end
        else
            evolip2=1;
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
            samp=moments{i}(c1,c2,t);
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
                dW{count}=sigma*randn(streams{count},[nnoise{count} 1]);
            end



            if ipdiag1 && ipdiag2
                cI1=evolip1.*c1;
                cI2=evolip2.*c2;
                ck1=evolip1.*G1(c1,c2,t,dW)*dt;%k1
                ck2=evolip2.*G2(c1,c2,t,dW)*dt;%k1
                c1=cI1+ck1/6;
                c2=cI2+ck2/6;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k2
                ck2=G2(ck1,ck2,t,dW)*dt;%k2
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k3
                ck2=G2(ck1,ck2,t,dW)*dt;%k3
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=evolip1.*(ck1+cI1);
                ck2=evolip2.*(ck2+cI2);
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k4
                ck2=G2(ck1,ck2,t,dW)*dt;%k4
                c1=evolip1.*c1+ck1/6;
                c2=evolip2.*c2+ck2/6;
            elseif ipdiag1
                cI1=evolip1.*c1;
                cI2=evolip2*c2;
                ck1=evolip1.*G1(c1,c2,t,dW)*dt;%k1
                ck2=evolip2*G2(c1,c2,t,dW)*dt;%k1
                c1=cI1+ck1/6;
                c2=cI2+ck2/6;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k2
                ck2=G2(ck1,ck2,t,dW)*dt;%k2
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k3
                ck2=G2(ck1,ck2,t,dW)*dt;%k3
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=evolip1.*(ck1+cI1);
                ck2=evolip2*(ck2+cI2);
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k4
                ck2=G2(ck1,ck2,t,dW)*dt;%k4
                c1=evolip1.*c1+ck1/6;
                c2=evolip2*c2+ck2/6;
            elseif ipdiag2
                cI1=evolip1*c1;
                cI2=evolip2.*c2;
                ck1=evolip1*G1(c1,c2,t,dW)*dt;%k1
                ck2=evolip2.*G2(c1,c2,t,dW)*dt;%k1
                c1=cI1+ck1/6;
                c2=cI2+ck2/6;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k2
                ck2=G2(ck1,ck2,t,dW)*dt;%k2
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k3
                ck2=G2(ck1,ck2,t,dW)*dt;%k3
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=evolip1*(ck1+cI1);
                ck2=evolip2.*(ck2+cI2);
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k4
                ck2=G2(ck1,ck2,t,dW)*dt;%k4
                c1=evolip1*c1+ck1/6;
                c2=evolip2.*c2+ck2/6;            
            else
                cI1=evolip1*c1;
                cI2=evolip2*c2;
                ck1=evolip1*G1(c1,c2,t,dW)*dt;%k1
                ck2=evolip2*G2(c1,c2,t,dW)*dt;%k1
                c1=cI1+ck1/6;
                c2=cI2+ck2/6;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k2
                ck2=G2(ck1,ck2,t,dW)*dt;%k2
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=ck1/2+cI1;
                ck2=ck2/2+cI2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k3
                ck2=G2(ck1,ck2,t,dW)*dt;%k3
                c1=c1+ck1/3;
                c2=c2+ck2/3;
                ck1=evolip1*(ck1+cI1);
                ck2=evolip2*(ck2+cI2);
                t=t+dt/2;
                ck1=G1(ck1,ck2,t,dW)*dt;%k4
                ck2=G2(ck1,ck2,t,dW)*dt;%k4
                c1=evolip1*c1+ck1/6;
                c2=evolip2*c2+ck2/6;
            end


            %do filters
            for i=1:length(filters)
                [c1,c2]=filters{i}(c1,c2,t);
            end


            stepcount=stepcount+1;
            takeSample=~mod(stepcount,sint*2);%half-step samples
            for i=1:length(takeSample)
                if takeSample(i)
                    samp=moments{i}(c1,c2,t);

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
        
    else
        sdata=sdataf;
        tdata=tdataf;
    end
    
end

