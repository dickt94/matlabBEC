load('test.mat')
for i=1:length(ks)
    i
    k=double(ks(i));
    testx=xs{i}.';
    testw=ws{i}.';
    [x,w]=gausshermite(k);
    norm(x-testx)
    norm(w-testw)
    if norm(x-testx)>k*eps('double')||norm(w-testw)>k*eps('double')
        error('output does not match test output for k=%i within precision', k);
    end
end