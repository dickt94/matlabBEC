function [samples,times]=spgpesample()
    %SPGPE simulation in HG basis, to produce initial states for integration of
    %NPW simulation.
    nmodes=128;
    c0=zeros([nmodes 1]);
    c0(1)=1;%start in ground state

    nx=(0:nmodes-1).';

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

    size(w_4f)

    function F=f(c,t)

        F=-(1i+gamma_x)*(-mu).*c;

        %add the nonlinear term
        psi=trans_4f*c;
        F=F-g_1D*(1i+gamma_x)*invtrans_4f*(w_4f.*conj(psi).*psi.*psi);

    end

    function G=g1(c,t)
        G=1/sqrt(2)*sgT*ones(nmodes,1);
    end

    function G=g2(c,t)
        G=1i/sqrt(2)*sgT*ones(nmodes,1);
    end

    function s=norm(c,t)
        s=c'*c;
    end

    function d=dens(c,t)
        d=conj(c).*c;
    end

    ipevol=-(1i+gamma_x)*(nx+0.5);

    [samples,times]=rk4int(c0,@f,ipevol,true,2,{@g1,@g2},[nmodes nmodes],[103029 837189039],0,time_int,50000,{@norm,@dens},[500 500],[true true],{});
end