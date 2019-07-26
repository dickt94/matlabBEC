function [samples,times]=spgpesample()
    %SPGPE simulation in HG basis, to produce initial states for integration of
    %NPW simulation.
    nmodes=20;
    c0=zeros([nmodes 1]);
    %desired particle number
    npart=1000;
    c0(1)=sqrt(npart)*1;%start in (HG) ground state
    
    %todo: start at TF ground state

    nx=(0:nmodes-1).';

    % Effective 1D interaction strength
    g_1D = 0.1;%2.0 * Lambda * a_s / sqrt( hbar / (m * omega_x) );

    % Chemical potential
    mu = 40;

    % EVOLUTION

    % Time interval of integration
    time_int = 150.0;


    %generate transforms here
    [x_4f,w_4f,trans_4f]=nfieldtrans(nmodes,4);
    invtrans_4f=trans_4f';

    function F=f(c,t)

        F=0;

        %add the nonlinear term
        psi=trans_4f*c;
        F=F-g_1D*invtrans_4f*(w_4f.*conj(psi).*psi.*psi);

    end

    function c=renorm(c,~)
        norm=c'*c;
        c=c/sqrt(norm)*sqrt(npart);
    end

    function s=norm(c,~)
        s=c'*c;
    end

    function d=dens(c,~)
        d=conj(c).*c;
    end

    function smp=fieldsamp(c,~)
        smp=c;
    end

    ipevol=-(nx+0.5);

    [samples,times]=rk4int(c0,@f,ipevol,true,0,{},{},[],0,time_int,30000000,{@norm,@dens,@fieldsamp},[500 500 10000],{@renorm},false);
end