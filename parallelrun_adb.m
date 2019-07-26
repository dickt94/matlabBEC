function parallelrun(npaths)
    rng('shuffle');
    parfor i=1:npaths
        [S,T]=npw_adiabaticramp;
        data{i}=S;
        times{i}=T;
    end
    
    save('initsim.mat','data','times','-v7.3');

    
end