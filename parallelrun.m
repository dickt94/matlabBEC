function parallelrun(npaths)
    parfor i=1:npaths
    
        [S,T]=npw_nofilter;
        data{i}=S;
        times=T;
    end
    
    save('output.mat','data');
    save('output.mat','times');
    
end