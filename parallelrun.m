function parallelrun(npaths)
    rng('shuffle');
    load('initsim.mat','data')
    parfor i=1:npaths
        [S,T]=npw_nofilter(data{i}{2}(:,:,end));
        data{i}=S;
        times{i}=T;
    end
    
    save('output.mat','data','times','-v7.3');

    
end