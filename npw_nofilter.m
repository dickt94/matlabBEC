function [samples,times]=npw_nofilter()
    %NPW simulation without system-filter separation
    nmodes=128;
    npaths=2;
    initw=zeros([1 npaths]);%initial log weights
    
    %import initial state
    %field representation should be c-matrix of size [npaths nmodes+1].
    %last element of each row holds the log weight.
    sample=zeros([nmodes npaths]);
    load('testsample.mat','sample');
    c0=[sample(:,1:npaths);initw];
    
    nx=(0:nmodes-1).';
    
    sqrtn12=sqrt((nx+1)/2);

    % FIXED CONSTANTS - ONLY NEEDED FOR EFFECTIVE INTERACTION ENERGY
    hbar = 1.05457173e-34;
    m = 86.909180527*1.66053892e-27;
    a_s = 5.313e-9;
    omega_x = 2.0 * pi * 20.0;
    % 1D BOSE GAS PARAMETERS (all in H.O. units)

    % Ratio of trapping frequencies (omega_perp / omega_x)
    Lambda = 50.0;

    % Temperature
    T = 50.0;

    % Effective 1D interaction strength
    g_1D = 2.0 * Lambda * a_s / sqrt( hbar / (m * omega_x) );

    % Chemical potential
    mu = 0.8* T;

    % Growth rate (assumed to be constant accros grid). For equilibrium
    % distribution this precise value doesn't matter.
    gamma_x = 3.0e-2;
    % Prefactor for noise correlator. 
    sgT = sqrt(2.0 * gamma_x * T);

    % EVOLUTION

    % Time interval of integration
    time_int = 100.0;


    %generate transforms here
    [x_4f,w_4f,trans_4f]=nfieldtrans(nmodes,4);
    invtrans_4f=trans_4f';
    
    dumw=zeros([1 npaths]);
    function F=f(c,t)

        F=0;

        %add the nonlinear term
        psi=trans_4f*c(1:end-1,:); %only use the field part, not the weights
        F=[F-g_1D*(1i+gamma_x)*invtrans_4f*(w_4f.*conj(psi).*psi.*psi);dumw];
        
        %feedback Hamiltonian
        
        %expect p
        cnp1=zeros(size(c(1:end-1,:)));
        cnp1(1:end-1,:)=c(2:end-1,:);
        
        size(sqrtn12)
        
        ep=2*sqrtn12.*(c(1:end-1,:)'*cnp1);
        

    end

    %real noise
    function G=g1(c,t,dW)
        G=[1/sqrt(2)*sgT*dW; 0];
    end

    function G=g2(c,t,dW)
        G=[1i/sqrt(2)*sgT*dW; 0];
    end

    %fictitious noises
    
    function s=norm(c,t)
        s=diag(c(1:end-1,:)'*c(1:end-1,:));
    end

    function d=dens(c,t)
        d=conj(c(1:end-1,:)).*c(1:end-1,:);
    end

    ipevol=[-(1i+gamma_x)*(nx+0.5-mu); 0];

    [samples,times]=rk4int(c0,@f,ipevol,true,2,{@g1,@g2},[nmodes nmodes],[103029 837189039],0,time_int,50000,{@norm,@dens},[500 500],{});
end