function singlerun(T,mu,alpha,res,npaths,outputDir,jobIndex)
    %ideally npaths should be number of cores on one NCI node
    rng('shuffle');
    clu=parcluster();
    clu.JobStorageLocation=[clu.JobStorageLocation num2str(jobIndex)];
    
    parpool(clu,npaths);
    parfor i=1:npaths
        [data{i},times{i}]=npw_adiabaticramp(T,mu);
        [data{i},times{i}]=npw_nofilter(data{i}{2}(:,:,end),alpha,res);
    end
    fname=[outputDir '/output-T' num2str(T) 'alpha' num2str(alpha) 'res' num2str(res) '.mat']
    
    save(fname,'data','times');
    
end
