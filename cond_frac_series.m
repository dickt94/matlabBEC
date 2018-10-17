function condfrac=cond_frac_series(fields,weights,ntimes,nmodes,npaths)
    condfrac=zeros([ntimes 1]);
    for count=1:ntimes
        psi=fields(:,:,count);
        w=exp(weights(:,count));
        w=w/sum(w);
        wpsi=sqrt(w).*psi;
        rho=wpsi'*wpsi-1/2*(nmodes/npaths)*diag(ones([npaths 1]));
        avg_N=sum(diag(rho));
        [V,D]=eig(rho);
        D2=sort(diag(D),'descend');
        condfrac(count)=D2(1)/avg_N;
    end
end