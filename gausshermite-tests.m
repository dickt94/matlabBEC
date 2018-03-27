load('test.mat')
ks
for i in range(length(ks))
    k=ks(i);
    testx=xs(i);
    testw=ws(i);
    [x,w]=gausshermite(k);
    x-testx
    w-testw
end