%tests a three-field integral that needs to be calculated for the NPW
%feedback simulation, with parameters that give a known analytic
%solution.

%this integral is taken from eq 74 of two_boson_spin_corr_notes.pdf by
%Stuart.
function res=integraltest()
    %generate a three-field grid with as many modes as possible - this
    %improves the approximation.
    
    [x,w,T]=nfieldtrans(371,3);
    
    hg1=zeros(371,1);
    hg2=zeros(371,1);
    m=0;n=0;
    hg1(n+1)=1;
    hg2(m+1)=1;
    phi1=T*hg1;
    phi2=T*hg2;
    
    psi=exp(-x.^2);
    
    plot(x,phi1,x,exp(-x.^2/3)*1/pi^(1/4));
    
    res=sum(w.*phi1.*phi2.*psi);
    
end