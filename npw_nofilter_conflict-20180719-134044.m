function [samples,times]=npw_nofilter()
    %NPW simulation without system-filter separation
    nmodes=50;
    npaths=100;
    initw=zeros([1 npaths]);%initial log weights
    
    %import initial state
    %field representation should be c-matrix of size [npaths nmodes+1].
    %last element of each row holds the log weight.
    sample=zeros([nmodes npaths]);
    load('testsample.mat','sample');
    c0=[sample(:,end-npaths+1:end);initw];

    
    %test the ground state
    %c0=zeros([nmodes+1 npaths]);
    %c0(1,:)=1;
    
    %test excited state
    %c0=zeros([nmodes+1 npaths]);
    %c0(1,end)=1;
    
    nx=(0:nmodes-1).';
    
    sqrtn12=sqrt((nx+1)/2);
    sqrtn2=sqrt(nx/2);
    
    sqrtnp1=sqrt(nx+1);
    sqrtnp2=sqrt(nx+2);
    
    sqrtn=sqrt(nx);
    sqrtnm1=sqrt(nx-1);

    % FIXED CONSTANTS - ONLY NEEDED FOR EFFECTIVE INTERACTION ENERGY
    hbar = 1.05457173e-34;
    m = 86.909180527*1.66053892e-27;
    a_s = 5.313e-9;
    omega_x = 2.0 * pi * 20.0;
    % 1D BOSE GAS PARAMETERS (all in H.O. units)

    % Ratio of trapping frequencies (omega_perp / omega_x)
    Lambda = 50.0;
    
    % Effective 1D interaction strength
    g_1D = 2.0 * Lambda * a_s / sqrt( hbar / (m * omega_x) );

    % Chemical potential
    %mu = 0;%0.8* T;

    % Growth rate (assumed to be constant accros grid). For equilibrium
    % distribution this precise value doesn't matter.
    %3.0e-2;
    % Prefactor for noise correlator. 

    
    %measurement params
    %strength
    alpha=0.001;
    sqrta=sqrt(alpha);
    %resolution
    r=0.01;
    %efficiency
    eta=1.0;
    sqrte=sqrt(eta);
    
    %feedback params
    fbsl=1.0;
    fbbr=0.25;
    fbnl=200.0;
    
    % EVOLUTION

    % Time interval of integration
    time_int = 100.0;


    %generate transforms here
    [x_4f,w_4f,trans_4f]=nfieldtrans(nmodes,4);
    invtrans_4f=trans_4f';
    [x_2f,w_2f,trans_2f]=nfieldtrans(nmodes,2);
    invtrans_2f=trans_2f';
    [x_3f,w_3f,trans_3f]=nfieldtrans(nmodes,3);
    invtrans_3f=trans_3f';
    %dummy initial state
    %psi0=exp(-(x_2f-1).^2 /2);
    
    %c0=[invtrans_2f*psi0;0]
    
    dumw=zeros([1 npaths]);
    %initnorm2=diag(c0(1:end-1,:)'*c0(1:end-1,:)); 
    %function c=renormalise(c,t)
    %    norm2=diag(c(1:end-1,:)'*c(1:end-1,:));
    %    for path=1:npaths
    %        c(1:end-1,path)=c(1:end-1,path)/sqrt(norm2(path)/initnorm2(path));
    %    end
    %end
        
    sqrt2p=sqrt(2*pi);

    %delta(x,x) for position density
    delt2f=trans_2f*trans_2f';
    
    %precompute Fourier transform of truncated delta - four-field
    
    Fdelt_4f=trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*sum(trans_3f.^2,2)));
    %define kernel in four-field grid
    nu_k=sqrt(r/(2*gamma(5/4)))*exp(-(r*x_4f/sqrt(2)).^4/2);%k scaled appropriately for the four-field grid
    
    zeta=nu_k.*Fdelt_4f;
    %precompute position integral of zeta^2
    int_zeta2=sum(w_4f.*zeta.*zeta);
    %zeta2=zeta.*zeta;
    
    function F=f(c,t)

        F=0;

        %add the nonlinear term
        psi=trans_4f*c(1:end-1,:); %only use the field part, not the weights
        F=F-[g_1D*(1i)*invtrans_4f*(w_4f.*conj(psi).*psi.*psi);dumw];
        
        %feedback Hamiltonian
        
      
        %cnp1=[c(2:end-1,:);zeros([1 npaths])];
        %cnm1=[zeros([1 npaths]); c(1:end-2,:)];
        %cnp2=[c(3:end-1,:);zeros([2 npaths])];
        %cnm2=[zeros([2 npaths]); c(1:end-3,:)];
        
        w=exp(c(end,:)).';
        w=w/sum(w);
        
        %norm
        %norm=0;
        %for path=1:npaths
        %    norm=norm+sum(w(path)*conj(c(1:end-1,path)).*c(1:end-1,path));
        %end
        norm=sum((conj(c(1:end-1,:)).*c(1:end-1,:)-0.5*ones([nmodes npaths]))*w);
        
        
        %expect p and xp+px
        %ep=0;
        %flat= @(M) M(:);
        ep=2/norm*sqrtn12.'*imag(conj(c(1:end-1,:)).*[c(2:end-1,:);zeros([1 npaths])])*w;
        %ep=2/norm*sum(flat(w.'.*conj(c(1:end-1,:)).*[c(2:end-1,:);zeros([1 npaths])].*sqrtn12));        
        xppx=2/norm*(sqrtnp1.*sqrtnp2).'*imag(conj(c(1:end-1,:)).*[c(3:end-1,:);zeros([2 npaths])])*w;
        
        %calc expected values with loops - obsolete and slower
        %xppx=0;
        %for path=1:npaths
            %ep=ep+2*w(path)*sum(imag(conj(c(1:end-1,path)).*[c(2:end-1,path);0].*sqrtn12))/norm;
        %    xppx=xppx+2*w(path)*sum(imag(conj(c(1:end-1,path)).*sqrtnp1.*sqrtnp2.*[c(3:end-1,path);0;0]))/norm;
        %end
        %ep=diag(2*imag(c(1:end-1,:)'*(sqrtn12.*[c(2:end-1,:);zeros([1 npaths])]))).'*w/norm;
        
        %expect xp+px
        
        %xppx=2*diag(imag(c(1:end-1,:)'*(sqrtnp1.*sqrtnp2.*[c(3:end-1,:);zeros([2 npaths])]))).'*w/norm;
        
        %add feedback
        F=F-1i*[fbsl*ep*(sqrtn2.*[zeros([1 npaths]); c(1:end-2,:)]+sqrtn12.*[c(2:end-1,:);zeros([1 npaths])])+fbbr*xppx*(0.5*(sqrtn.*sqrtnm1.*[zeros([2 npaths]); c(1:end-3,:)]+sqrtnp1.*sqrtnp2.*[c(3:end-1,:);zeros([2 npaths])])+(nx+0.5).*c(1:end-1,:));dumw];
        
        %energy-damping feedback
        F=F+1i*[fbnl/norm*invtrans_4f*(w_4f.*psi.*imag(sum(w.'.*(conj(psi).*(trans_4f*(0.5*(sqrtnp1.*sqrtnp2.*[c(3:end-1,:);zeros([2 npaths])]+sqrtn.*sqrtnm1.*[zeros([2 npaths]);c(1:end-3,:)])-(nx+0.5).*c(1:end-1,:)))),2)));dumw];
        
        %calculate deterministic weight bits - notation from my notes
        %2018/02/01
        psi_3f=trans_3f*c(1:end-1,:);
        xi_j=nu_k.*trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psi_3f).*psi_3f));
        
        %for loop is probably slow - replace with matrix algebra if
        %possible
%         m_pj=zeros([npaths npaths]);
%         for countp=1:npaths
%            for countj=1:npaths
%                xi_jj=flip(xi_j(:,countj));
%                xi_jp=xi_j(:,countp);
%                m_pj(countp,countj)=2*pi*sum(w_4f.*(xi_jp.*xi_jj-zeta.*(xi_jj+xi_jp)+zeta2));
%            end
%         end
%         
%         m_pjb=m_pj;
        
        %compute m_pj using matrices
        m_pj=2*pi*(xi_j.'*(w_4f.*flip(xi_j,1))+int_zeta2-xi_j.'*(w_4f.*zeta)-(w_4f.*zeta).'*flip(xi_j,1));
        
        %m_j^2 is the diagonal of m_pj
        %calculate the deterministic weight evolution
        wmmj=sum(w.*m_pj,1);
        wmwm=sum(w.'.*wmmj);
        mj2=diag(m_pj).';
        wmj2=sum(w.'.*mj2);
        detw=2*eta*alpha*(2*(wmmj-wmwm)-(mj2-wmj2));
        F=F+[zeros([nmodes npaths]);detw];
        
    end

    %fictitious noises
    %calculate functions xi_n(x)
    phi_nk=((-1i).^(nx.')).*trans_4f;
    d_mn=phi_nk'*(nu_k.*phi_nk);
    xi_nx=d_mn.'*invtrans_3f;

    function fnoise=fieldnoise(c,t,dW)
        %dW should be vector of [nmodes npaths] size, corresponding to
        %fictitious
        %noises in HG space
        fnoise=-1i*[sqrta*invtrans_3f*(w_3f.*((xi_nx.'*dW).*(trans_3f*c(1:end-1,:))));dumw];
    end
    %real noise
    
    %precompute Fourier transform of truncated delta - three-field
    
    Fdelt=trans_3f*(((-1i).^nx).*invtrans_3f*(w_3f.*sum(trans_3f.^2,2)));
    %define kernel with k scaled to 3-field grid
    nu_k3f=sqrt(r/(2*gamma(5/4)))*exp(-(r*x_3f/sqrt(3/2)).^4/2);%k scaled appropriately for the three-field grid
    
    function wev=weightevol(c,t,rn)
        %rn is vector of [nmodes 1] size, corresponding to real noise in HG
        %space
        psi=trans_3f*c(1:end-1,:);
        %project Fourier transform of |psi|^2 into HG basis using three-field
        Fpsi2=trans_3f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psi).*psi));
        %fourier noise
        knoise=trans_3f*((1i.^nx).*rn);
        %integrate to get the evolution term
        noiseterm_j=2*sqrte*sqrta*sqrt2p*sum(w_3f.*nu_k3f.*knoise.*(Fpsi2-Fdelt),1);
        w=exp(c(end,:));
        w=w/sum(w);
        noiseterm_j=noiseterm_j-sum(w.*noiseterm_j);
        wev=[zeros([nmodes npaths]);noiseterm_j];
    end
    
    %breeding filter
    breedgap=12;
    breedcount=0;
    function c=breed(c,t)
        
        [maxw,maxi]=max(c(end,:));
        [minw,mini]=min(c(end,:));
        
        %breed
        while maxw-minw > breedgap
            fprintf("breeding\n");
            c(end,mini)=maxw-log(2);
            c(end,maxi)=maxw-log(2);
            c(1:end-1,mini)=c(1:end-1,maxi);
            [maxw,maxi]=max(c(end,:));
            [minw,mini]=min(c(end,:));
            breedcount=breedcount+1;
        end
        
        %normalise
        c(end,:)=c(end,:)-maxw;
        
    end


    function s=norm(c,t)
        w=exp(c(end,:));
        w=w/sum(w);
        s=sum((conj(c(1:end-1,:)).*c(1:end-1,:)-0.5*ones([nmodes npaths]))*w);
    end

    function f=field(c,t)
        f=c;
    end

    delt_4f=
    function e=eE(c,t)
        %expected value of energy in HO units
        
        %linear
        w=exp(c(end,:));
        w=w/sum(w);
        dens_nx=conj(c(1:end-1,:)).*c(1:end-1,:)-1/2*ones([nmodes npaths]);
        e_lin=(nx+0.5).'*dens_nx*w.';
        
        %nl energy
        psi_4f=trans_4f*c(1:end-1,:); %only use the field part, not the weights
        %dens_x=conj(psi).*psi-1/2*diag(delt2f);
        enlpath=0.5*g_1D*sum(w_4f.*conj(psi).*conj(psi).*psi.*psi,1);
        e=e+w*enlpath'
        
    end

    function rho=obdm(c,t)
        w=exp(c(end,:));
        sqrtw=sqrt(w/sum(w));
        wpsi=sqrtw.*c(1:end-1,:);
        rho=wpsi*wpsi'-1/2*diag(ones(nmodes,1));
    end


    function rhoxy=obdmx(c,t)
        w=exp(c(end,:));
        sqrtw=sqrt(w/sum(w));
        psi=trans_2f*c(1:end-1,:);
        wpsi=sqrtw.*psi;
        rhoxy=wpsi*wpsi'-1/2*delt2f;
    end

    function e=evenness(c,t)
        w=exp(c(end,:));
        w=w/sum(w);
        e=1/sum(w.*w);
    end

    ipevol=[-(1i)*(nx+0.5); 0];

    [samples,times]=rk4int(c0,@f,ipevol,true,2,{@fieldnoise @weightevol},{[nmodes npaths],[nmodes 1]},[103029 837189039],0,time_int,50000,{@norm,@field,@eE,@obdm,@obdmx,@evenness},[500 100 500 100 100 500],{@breed},false);
    fprintf("bred %d times\n",breedcount);
end