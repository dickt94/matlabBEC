function [samples,times]=npw_gpuhf(initstate1,initstate2,alpha,measres)
    %initialise RNG (comment out for parfor, deal with random substreams
    %independently or use default behaviour
    %rng('shuffle');
    %NPW simulation without system-filter separation
    nmodes=size(initstate1,1)-1;
    npaths=size(initstate1,2);
    initw=zeros([1 npaths]);%initial log weights
    
    %import initial state
    %field representation should be c-matrix of size [npaths nmodes+1].
    %last element of each row holds the log weight.
    %sample=0;
    %load('sample-jl.mat','sample');
    %shuffle the sample
    %sample=sample(:,randperm(size(sample,2)));
    %take a path from the thermal SPGPE state and use it to generate an NPW
    %state
    
    %alpha0_n=sample(:,1);%[0;sqrt(1000);zeros(nmodes-2,1)];
%     n0k=zeros([nmodes npaths]);
%     phi0k=zeros([nmodes npaths]);
%     
%     for k=1:npaths
%        n0k(:,k)=poissrnd(abs(alpha0_n).^2);
%        for nlev=1:nmodes
%            if n0k(nlev,k)==0
%                phi0k(nlev,k)=2*pi*rand();
%            else
%                phi0k(nlev,k)=normrnd(angle(alpha0_n(nlev)),1/4*psi(1,n0k(nlev,k)+1));
%            end
%        end
%     end
    
    c01=initstate1;
    c02=initstate2;
    %c0(1:end-1,:)=repmat(alpha0_n,[1 npaths])+(randn(nmodes,npaths)+1i*randn(nmodes,npaths))/2;
    %c0(1:end-1,:)=sqrt(n0k+1/2).*exp(1i*phi0k);

    
    %test the ground state
    %c0=zeros([nmodes+1 npaths]);
    %c0(1,:)=1;
    
    %test excited state
    %c0=zeros([nmodes+1 npaths]);
    %c0(1,end)=1;
    
    
    nx=gpuArray(complex((0:nmodes-1).'));
    
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
    g_1D = 0.01;%0.1;%0.4;%2.0 * Lambda * a_s / sqrt( hbar / (m * omega_x) );

    % Chemical potential
    %mu = 0;%0.8* T;

    % Growth rate (assumed to be constant accros grid). For equilibrium
    % distribution this precise value doesn't matter.
    %3.0e-2;
    % Prefactor for noise correlator. 

    
    %measurement params
    %strength
    sqrta=sqrt(alpha);
    %resolution
    r=measres;
    %efficiency
    eta=1;
    sqrte=sqrt(eta);
    
    %feedback params
    fbsl=1.0;%1.0;
    fbbr=0;%0.25;
    fbnl=0.01;%100.0;
    
    % EVOLUTION

    % Time interval of integration
    time_int=400.0;
    nsteps=30000;


    %generate transforms here
    [x_4f,w_4f,trans_4f]=nfieldtrans(nmodes,4);
    x_4f=gpuArray(x_4f);
    w_4f=gpuArray(w_4f);
    trans_4f=gpuArray(trans_4f);
    invtrans_4f=trans_4f';
    [x_2f,w_2f,trans_2f]=nfieldtrans(nmodes,2);
    x_2f=gpuArray(x_2f);
    w_2f=gpuArray(w_2f);
    trans_2f=gpuArray(trans_2f);
    invtrans_2f=trans_2f';
    [x_3f,w_3f,trans_3f]=nfieldtrans(nmodes,3);
    x_3f=gpuArray(x_3f);
    w_3f=gpuArray(w_3f);
    trans_3f=gpuArray(trans_3f);
    invtrans_3f=trans_3f';
    
    %three-field fourier transform - useful for things like FT of density
    f3f=trans_3f*(((-1i).^nx).*invtrans_3f);
    fourier_3f=w_3f'.*f3f;
    invfourier_3f=w_3f'.*f3f';
    
    %dummy initial state
    %psi0=exp(-(x_2f-1).^2 /2);
    
    %c0=[invtrans_2f*psi0;0]
    
    %kick the initial state
    %psi0=trans_2f*c0(1:end-1,:);
    %k0=0;
    %psi0=exp(1i*k0*x_2f).*psi0;
    %c0(1:end-1,:)=invtrans_2f*psi0;
    
    
    %sample a coherent BEC (for testing against Michael's NJP)
    %coherent state magnitudes alpha
    %part_num=1124;
    %x0=0;
    %sig=4.0;
    %alpha0_x=sqrt(part_num/(sqrt(2*pi)*sig))*exp(-(x_2f-x0).^2/(4*sig^2));
    %transform the coherent state magnitudes to HG space
    %alpha0_n=invtrans_2f*(w_2f.*alpha0_x);
    %n0k=zeros([nmodes npaths]);
    %phi0k=zeros([nmodes npaths]);
    
    %for k=1:npaths
    %    n0k(:,k)=poissrnd(abs(alpha0_n).^2);
    %    for nlev=1:nmodes
    %        if n0k(nlev,k)==0
    %            phi0k(nlev,k)=2*pi*rand();
    %        else
    %            phi0k(nlev,k)=normrnd(angle(alpha0_n(nlev)),1/4*psi(1,n0k(nlev,k)+1));
    %        end
    %    end
    %end
    
    %c0=zeros([nmodes+1 npaths]);
    %c0(1:end-1,:)=sqrt(n0k+1/2).*exp(1i*phi0k);
    
    %super dumb initial state
    %c0(1:end-1)=alpha0_n;
    
    %'ground' state
    %load('gssample.mat','alpha0_n');
    %alpha0_n=zeros(nmodes,1);
    %alpha0_n(1)=sqrt(1000);
    %c0(1:end-1,:)=repmat(alpha0_n,[1 npaths])+(randn(nmodes,npaths)+1i*randn(nmodes,npaths))/2;

    
    dumw=gpuArray(zeros([1 npaths]));
    %initnorm2=diag(c0(1:end-1,:)'*c0(1:end-1,:)); 
    %function c=renormalise(c,~)
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
    nu_k=sqrt(r)/(2*gamma(5/4))*exp(-(r*x_4f/sqrt(2)).^4/2);%k scaled appropriately for the four-field grid
    
    %define kernel with k scaled to 3-field grid
    nu_k3f=sqrt(r)/(2*gamma(5/4))*exp(-(r*x_3f/sqrt(3/2)).^4/2);%k scaled appropriately for the three-field grid
    
    nu_k2=nu_k.^2;
    nu_k3f2=nu_k3f.^2;
    
    zeta=nu_k.*Fdelt_4f;
    %precompute position integral of zeta^2
    int_zeta2=sum(w_4f.*zeta.*zeta);
    %zeta2=zeta.*zeta;
    
    delt_4f=diag(trans_4f*trans_4f');
    
    
    function F=ffil(csys,c,~)
        F=0;
        %add the nonlinear term
        psi=trans_4f*c;
        norm=c'*c;
        
        F=F-(g_1D)*(1i)*invtrans_4f*(w_4f.*(abs(psi).^2.*psi));
        
        %feedback Hamiltonian
        ep=2/norm*sum(sqrtn12.'*imag(conj(c).*[c(2:end);0]));
        
        xppx=2/norm*sum((sqrtnp1.*sqrtnp2).'*imag(conj(c).*[c(3:end);0;0]));
        F=F-1i*fbsl*ep*(sqrtn2.*[0; c(1:end-1)]+sqrtn12.*[c(2:end);0])+fbbr*xppx*(0.5*(sqrtn.*sqrtnm1.*[0;0;c(1:end-2)]+sqrtnp1.*sqrtnp2.*[c(3:end);0;0])+(nx+0.5).*c);
        
        
        %energy-damping control
        ved=imag(conj(psi).*(trans_4f*(0.5*(sqrtnp1.*sqrtnp2.*[c(3:end);0;0]+sqrtn.*sqrtnm1.*[0;0;c(1:end-2)])-(nx+0.5).*c)));
        F=F+1i*fbnl/norm*invtrans_4f*(w_4f.*ved.*psi);
        
        psi_3f=trans_3f*c;
        %F=F+4*eta*alpha*invtrans_3f*(w_3f.*psi_3f.*trans_3f*((-1i).^nx.*invtrans_3f*(w_3f.*nu_k3f.^2.*trans_3f*((1i).^nx.*invtrans_3f*(w_3f.*abs(psi_3f).^2)))));
        
        %deterministic part of innovations
        ws=exp(csys(end,:)).';
        ws=ws/sum(ws);
        psis_3f=trans_3f*csys(1:end-1,:);
        denss=(conj(psis_3f).*psis_3f)*ws;
        F=F+4*eta*alpha*invtrans_3f*(w_3f.*psi_3f.*trans_3f*((-1i).^nx.*invtrans_4f*(w_4f.*nu_k2.*trans_4f*((1i).^nx.*invtrans_3f*(w_3f.*denss)))));
        
    end

    function F=fsys(csys,c,~)

        F=0;
        
        %add the nonlinear term
        psisys=trans_4f*csys(1:end-1,:); %only use the field part, not the weights
        psi=trans_4f*c;
        
        w=exp(csys(end,:)).';
        w=w/sum(w);
        
        norm=c'*c;
        %nl feedback
        %ufb=-sum(w_4f.*imag(sum(w.'.*conj(psi).*abs(psi).^2.*(trans_4f*(0.5*(sqrtnp1.*sqrtnp2.*[c(3:end-1,:);zeros([2 npaths])]+...
        %    sqrtn.*sqrtnm1.*[zeros([2 npaths]);c(1:end-3,:)])-(nx+0.5).*c(1:end-1,:))),2)))/(norm^2);
        %wigner correction added 30/08
        
        
        F=F-[(g_1D)*(1i)*invtrans_4f*(w_4f.*(abs(psisys).^2.*psisys-delt_4f.*psisys));dumw];%-delt_4f.*psi
        
        %feedback Hamiltonian
        
      
        %cnp1=[c(2:end-1,:);zeros([1 npaths])];
        %cnm1=[zeros([1 npaths]); c(1:end-2,:)];
        %cnp2=[c(3:end-1,:);zeros([2 npaths])];
        %cnm2=[zeros([2 npaths]); c(1:end-3,:)];
        
        %norm
        %norm=0;
        %for path=1:npaths
        %    norm=norm+sum(w(path)*conj(c(1:end-1,path)).*c(1:end-1,path));
        %end
        
        
        %expect p and xp+px, for filter!
        %ep=0;
        %flat= @(M) M(:);
        ep=2/norm*sum(sqrtn12.'*imag(conj(c).*[c(2:end);0]));
        %ep=2/norm*sum(flat(w.'.*conj(c(1:end-1,:)).*[c(2:end-1,:);zeros([1 npaths])].*sqrtn12));        
        xppx=2/norm*sum((sqrtnp1.*sqrtnp2).'*imag(conj(c).*[c(3:end);0;0]));
        
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
        F=F-1i*[fbsl*ep*(sqrtn2.*[zeros([1 npaths],'gpuArray'); csys(1:end-2,:)]+sqrtn12.*[csys(2:end-1,:);zeros([1 npaths],'gpuArray')])+fbbr*xppx*(0.5*(sqrtn.*sqrtnm1.*[zeros([2 npaths],'gpuArray'); csys(1:end-3,:)]+sqrtnp1.*sqrtnp2.*[csys(3:end-1,:);zeros([2 npaths],'gpuArray')])+(nx+0.5).*csys(1:end-1,:));dumw];
        
        %energy-damping feedback
        %compute energy-damping feedback potential
        ved=imag(conj(psi).*(trans_4f*(0.5*(sqrtnp1.*sqrtnp2.*[c(3:end);0;0]+sqrtn.*sqrtnm1.*[0;0;c(1:end-2)])-(nx+0.5).*c)));
        F=F+1i*[fbnl/norm*invtrans_4f*(w_4f.*ved.*psisys);dumw];
        
        
        %compute 'quantum noise control' (from michael's NJP 2013)
        psisys_3f=trans_3f*csys(1:end-1,:);
        %dpsi_3f_dx=trans_3f*(sqrtn12.*[c(2:end-1,:);zeros([1 npaths])]+sqrtn2.*[zeros([1 npaths]);c(1:end-2,:)]);
        %compute potential using fourier transforms
        %vqnc=invfourier_3f*(1i*x_3f.*(fourier_3f*imag(psi_3f.*dpsi_3f_dx)).*nu_k3f2)*w;
        %F=F-1i*fbnl/norm*[invtrans_3f*(vqnc.*psi_3f);dumw];
        
        %calculate deterministic weight bits - notation from my notes
        %2018/02/01
        xi_j=nu_k.*trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psisys_3f).*psisys_3f));
        
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
        %
        
        %m_j^2 is the diagonal of m_pj
        %calculate the deterministic weight evolution
        wmmj=sum(w.*m_pj,1);
        wmwm=sum(w.'.*wmmj);
        mj2=diag(m_pj).';
        wmj2=sum(w.'.*mj2);
        detw=+2*eta*alpha*(2*(wmmj-wmwm)-(mj2-wmj2));
        
        F=F+[zeros([nmodes npaths],'gpuArray');detw];
        
    end



    %fictitious noises
    %calculate functions xi_n(x)
    phi_nk=((-1i).^(nx.')).*trans_4f;
    d_mn=phi_nk'*(nu_k.*phi_nk);
    xi_nx=d_mn.'*invtrans_3f;

    function fnoise=fieldnoise(c,dW)
        %dW should be vector of [nmodes npaths] size, corresponding to
        %fictitious
        %noises in HG space
        fnoise=-1i*[sqrta*invtrans_3f*(w_3f.*((xi_nx.'*dW).*(trans_3f*c(1:end-1,:))));dumw];
    end
    
    %fictitious noises act only on the system if the filter is HF
    fieldnoise1=@(c1,~,~,dW) fieldnoise(c1,dW);


    %real noise - the stochastic evolution is the same for system and
    %filter
    
    %precompute Fourier transform of truncated delta - three-field
    
    Fdelt=trans_3f*(((-1i).^nx).*invtrans_3f*(w_3f.*sum(trans_3f.^2,2)));
    
    function wev=weightevol(c,rn)
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

    weightevol1=@(c1,~,~,dW) weightevol(c1,dW);
    
    
    function wev=weightevol2(~,c,~,rn)
        %innovations noise for the filter
        psi=trans_3f*c;
        knoise=trans_3f*((1i.^nx).*rn);
        
        wev=2*sqrte*sqrta*invtrans_3f*(w_3f.*psi.*trans_3f*(((-1i).^nx).*invtrans_3f*(w_3f.*nu_k3f.*knoise)));
    end
    
    %breeding filter
    breedgap=10;
    breedcount=0;
    function [c1,c2]=breed(c1,c2,~)
        
        [maxw,maxi]=max(c1(end,:));
        [minw,mini]=min(c1(end,:));
        
        %breed
        while maxw-minw > breedgap
            c1(end,mini)=maxw-log(2);
            c1(end,maxi)=maxw-log(2);
            c1(1:end-1,mini)=c1(1:end-1,maxi);
            [maxw,maxi]=max(c1(end,:));
            [minw,mini]=min(c1(end,:));
            breedcount=breedcount+1;
        end
        
        %normalise
        c1(end,:)=real(c1(end,:)-maxw);
        
        %normalise HF filter
        c2=c2/sqrt(c2'*c2)*sqrt(2000);
        
    end


    %todo adapt sampling for system-filter separation. Get rid of debug
    %sampling for the whole complex field with noises - moments only.

    function s=norm(c)
        w=exp(c(end,:));
        w=w/sum(w);
        s=sum((conj(c(1:end-1,:)).*c(1:end-1,:)-0.5*ones([nmodes npaths]))*w.');
    end

    normS=@(c,~,~) norm(c);
    normF=@(~,c,~) c'*c;

    delt4f=trans_4f*trans_4f';
    function e=eE(c)
        %expected value of energy in HO units
        
        %linear
        w=exp(c(end,:));
        w=w/sum(w);
        dens_nx=conj(c(1:end-1,:)).*c(1:end-1,:)-1/2*ones([nmodes npaths]);
        e=(nx+0.5).'*dens_nx*w.';
        
        %nl energy
        psi=trans_4f*c(1:end-1,:); %only use the field part, not the weights
        
        %dens_x=conj(psi).*psi-1/2*diag(delt4f);
        %enlpath=0.5*g_1D*sum(w_4f.*(dens_x.^2-dens_x.*diag(delt4f)),1);
        %proper wigner correction (29/8/2018)
        dens_x_sym=conj(psi).*psi;
        enlpath=0.5*g_1D*sum(w_4f.*(dens_x_sym.^2-2*dens_x_sym.*diag(delt4f)+0.5*diag(delt4f).^2),1);
        e=e+w*enlpath';
        %energy per particle
        avg_N=sum((conj(c(1:end-1,:)).*c(1:end-1,:)-0.5*ones([nmodes npaths]))*w.');
        e=e/avg_N;
        
    end
    eS=@(c,~,~) eE(c);
    
    %hartree-fock energy
    function e=eEF(c)
        e=(nx+0.5).'*(conj(c).*c);
        psi=trans_4f*c;
        e=e+0.5*g_1D*w_4f.'*(abs(psi).^4);
        
        e=e/(c'*c);
    end
    
    eF=@(~,c,~) eEF(c);

    function rho=obdm(c)
        w=exp(c(end,:));
        sqrtw=sqrt(w/sum(w));
        wpsi=sqrtw.*c(1:end-1,:);
        rho=wpsi*wpsi'-1/2*diag(ones(nmodes,1));
    end

    rhoS=@(c,~,~) obdm(c);
    
    rhoF=@(~,c,~) c*c';

    function rhoxy=obdmx(c)
        w=exp(c(end,:));
        sqrtw=sqrt(w/sum(w));
        psi=trans_2f*c(1:end-1,:);
        wpsi=sqrtw.*psi;
        rhoxy=wpsi*wpsi'-1/2*delt2f;
    end

    rhoxS=@(c,~,~) obdmx(c);
    rhoxF=@(~,c,~) trans_2f*c*(trans_2f*c)';

%    transk_2f=((-1i).^nx).'.*trans_2f;
%    function rhokxky=obdmk(c,~)
%        w=exp(c(end,:));
%        sqrtw=sqrt(w/sum(w));
%        psi=transk_2f*c(1:end-1,:);
%        wpsi=sqrtw.*psi;
%        rhokxky=wpsi*wpsi'-1/2*delt2f;
%    end

    function e=evenness(c)
        w=exp(c(end,:));
        w=w/sum(w);
        e=1/sum(w.*w);
    end
    evenS=@(c,~,~) evenness(c);
    
    function v=momentvar(c)
        w=exp(c(end,:)).';
        w=w/sum(w);
        %calculate the momentum variance
        p2=((nx+0.5).'*(conj(c(1:end-1,:)).*c(1:end-1,:))-(sqrtnp1.*sqrtnp2).'*real(conj(c(1:end-1,:)).*[c(3:end-1,:);zeros([2 npaths])]))*w;
        p=2*sqrtn12.'*imag(conj(c(1:end-1,:)).*[c(2:end-1,:);zeros([1 npaths])])*w;
        s=sum((conj(c(1:end-1,:)).*c(1:end-1,:)-0.5*ones([nmodes npaths]))*w);
        %wigner correction for p2
        p2=p2-sum(nx)/2;
        v=((p2/s)-(p/s)^2);
    end

    mvarS=@(c,~,~) momentvar(c);
    
    
    %todo test feedback moment convergence
    
    
    %calculate innovation signal - useful for showing
    %convergence
    
%     function sig=innovsig(csys,c,~)
%         psi_3f=trans_3f*c(1:end-1,:);
%         psis_3f=trans_3f*csys(1:end-1,:);
%         w=exp(c(end,:)).';
%         w=w/sum(w);
%         wsys=exp(csys(end,:)).';
%         wsys=wsys/sum(wsys);
%         
%         
%         xi_j=nu_k.*trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psi_3f).*psi_3f));
%         xi_sj=nu_k.*trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psis_3f).*psis_3f));
%         
%         %compute m_pj using matrices
%         m_pj=2*pi*(xi_j.'*(w_4f.*flip(xi_j,1))+int_zeta2-xi_j.'*(w_4f.*zeta)-(w_4f.*zeta).'*flip(xi_j,1));
%         
%         %compute matrix for sys-fil signal coupling
%         m_spj=2*pi*(xi_sj.'*(w_4f.*flip(xi_j,1))+int_zeta2-xi_sj.'*(w_4f.*zeta)-(w_4f.*zeta).'*flip(xi_j,1));
%         
%         %m_j^2 is the diagonal of m_pj
%         %calculate the deterministic weight evolution
%         wmmj=sum(w.*m_pj,1);
%         wmwm=sum(w.'.*wmmj);
%         
%         %add the deterministic bit of the measurement signal
%         wmsmj=sum(wsys.*m_spj,1);
%         wmswm=sum(wsys.'.*wmsmj);
%         sig=4*sqrte*(wmsmj-wmswm-wmmj+wmwm);
% 
%     end
    fieldsamp=@(c) c;
    fsS=@(c1,~,~) fieldsamp(c1);
    fsF=@(~,c2,~) fieldsamp(c2);
    
    function sig=meassig(c1,c2,~)
        
        w1=exp(c1(end,:)).';
        w1=w1/sum(w1);
        
        psis_3f=trans_3f*c1(1:end-1,:);
        psi_3f=trans_3f*c2;
        
        xi_j=trans_2f*(((1i).^nx).*invtrans_4f*(w_4f.*nu_k.*trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psi_3f).*psi_3f))));
        xi_sj=trans_2f*(((1i).^nx).*invtrans_4f*(w_4f.*nu_k.*trans_4f*(((-1i).^nx).*invtrans_3f*(w_3f.*conj(psis_3f).*psis_3f))));

        sig=(xi_sj*w1-xi_j);
         
    end

    function fm=filtermoments(~,c2,~)
        norm=c2'*c2;
        ex=2/norm*sqrtn12.'*real(conj(c2).*[c2(2:end);0]);
        ep=2/norm*sqrtn12.'*imag(conj(c2).*[c2(2:end);0]);
        xppx=2/norm*(sqrtnp1.*sqrtnp2).'*imag(conj(c2).*[c2(3:end);0;0]);
        
        fm=[ep;xppx;ex];
    end

    function sm=sysmoments(c1,~,~)
        w1=exp(c1(end,:).');
        w1=w1/sum(w1);
        ex=2/norm(c1)*sqrtn12.'*real(conj(c1(1:end-1,:)).*[c1(2:end-1,:);zeros([1 npaths])])*w1;
        ep=2/norm(c1)*sqrtn12.'*imag(conj(c1(1:end-1,:)).*[c1(2:end-1,:);zeros([1 npaths])])*w1;
        xppx=2/norm(c1)*(sqrtnp1.*sqrtnp2).'*imag(conj(c1(1:end-1,:)).*[c1(3:end-1,:);zeros([2 npaths])])*w1;
        ex2=1/norm(c1)*((nx+0.5).'*(conj(c1(1:end-1,:)).*c1(1:end-1,:))+2*(sqrtnp1.*sqrtnp2).'*real(c1(1:end-1,:).*[c1(3:end-1,:);zeros([2 npaths])]))*w1;
        sm=[ep;xppx;ex;ex2];
    end

    function fbmomentdiff=fbdiff(c1,c2,~)
        
        w1=exp(c1(end,:).');
        w1=w1/sum(w1);
        
        Nf=c2'*c2;
        
        ep=2/norm(c1)*sqrtn12.'*imag(conj(c1(1:end-1,:)).*[c1(2:end-1,:);zeros([1 npaths])])*w1-2/Nf*sqrtn12.'*imag(conj(c2).*[c2(2:end);0]);
        xppx=2/norm(c1)*(sqrtnp1.*sqrtnp2).'*imag(conj(c1(1:end-1,:)).*[c1(3:end-1,:);zeros([2 npaths])])*w1-2/Nf*(sqrtnp1.*sqrtnp2).'*imag(conj(c2).*[c2(3:end);0;0]);
        ex2=1/norm(c1)*((nx+0.5).'*(conj(c1(1:end-1,:)).*c1(1:end-1,:))+2*(sqrtnp1.*sqrtnp2).'*real(c1(1:end-1,:).*[c1(3:end-1,:);zeros([2 npaths])]))*w1...
            -1/Nf*((nx+0.5).'*(conj(c2).*c2)+2*(sqrtnp1.*sqrtnp2).'*(conj(c2).*[c2(3:end);0;0]));
        
        ex=2/norm(c1)*sqrtn12.'*real(conj(c1(1:end-1,:)).*[c1(2:end-1,:);zeros([1 npaths])])*w1-2/Nf*sqrtn12.'*real(conj(c2).*[c2(2:end);0]);
        
        
        fbmomentdiff=[ep;xppx;ex;ex2];
    end

    %find the correction to the filter from the innovations term
    function corr=filtercorr(c1,c2,~)
        
        w1=exp(c1(end,:)).';
        w1=w1/sum(w1);
        
        psi_3f=trans_3f*c2;
        psis_3f=trans_3f*c1(1:end-1,:);
        
        
        denss=(conj(psis_3f).*psis_3f)*w1;
        
        corr=2*eta*alpha*invtrans_3f*(w_3f.*psi_3f.*trans_3f*((1i).^nx.*invtrans_4f*(w_4f.*nu_k2.*trans_4f*((-1i).^nx.*invtrans_3f*(w_3f.*denss)))));
        
         
    end

    function diff=densdiff(cs,cf,~)
        ws=exp(cs(end,:)).';
        ws=ws/sum(ws);
        psif_2f=trans_2f*cf;
        psis_2f=trans_2f*cs(1:end-1,:);
        diff=sum(w_2f.*abs(conj(psif_2f).*psif_2f-(conj(psis_2f).*psis_2f)*ws));
    end

%    function noisesamp=ns(dW)
%        noisesamp=dW;
%    end

%    function ved=fbV(c,~)
%        psi=trans_4f*c(1:end-1,:);
%        w=exp(c(end,:));
%        w=w/sum(w);
%        ved=imag(sum(w.*conj(psi).*(trans_4f*(0.5*(sqrtnp1.*sqrtnp2.*[c(3:end-1,:);zeros([2 npaths])]+sqrtn.*sqrtnm1.*[zeros([2 npaths]);c(1:end-3,:)])-(nx+0.5).*c(1:end-1,:))),2));
%    end

    ipevol=[-(1i)*(nx+0.5); 0];
    ipevol2=-(1i)*(nx+0.5);

    seed1=mod(randi([0 nsteps*nmodes*npaths]),2^32);
    seed2=mod(randi([0 nsteps*nmodes*npaths]),2^32);
    
    dumfunc=@(~,~,~,~) 0;
    
    [samples,times]=rk4int_gpu(c01,c02,@fsys,@ffil,ipevol,ipevol2,true,true,2,{fieldnoise1 weightevol1},{dumfunc @weightevol2},{[nmodes npaths],[nmodes 1]},[seed1 seed2],0,time_int,nsteps,...
        {normS,normF,eS,eF,rhoS,rhoF,rhoxS,rhoxF,evenS,fsS,fsF,@meassig,@fbdiff,@filtermoments,@sysmoments,@filtercorr,@densdiff},...
        [500 500 500 500 100 100 100 100 500 1 1 100 500 500 500 100 500],{@breed},false);
    fprintf("bred %d times\n",breedcount);
end