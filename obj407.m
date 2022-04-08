function object=obj407(dict,e,sig)
    espan = {[0 1] [0 1] [0 1] [0 1] [0 1] [0 1] [0 25] [0 25] [0 25] [0 25]}; 
    T=[198 248 273 298 348 398 298 298 298 298];
    srate=[3e-05 3e-05 3e-05 3e-05 3e-05 3e-05 4e-05 3e-04 3e-03 2e-02];
    sig_ss=[1.9 0.99 0.71 0.6 0.49 0.37 0.48 0.62 0.97 1.44];
    chi_0=0.25;
    G_298 = 2.83e3; 
    c0 = dict(1); 
    c1 = dict(2); 
    Kx = c0*exp(T/c1);
    K_pf=dict(3);
    t_0=dict(4);
    T_p=dict(5);
    sig_T = zeros(size(sig_ss)); 
    mu_T = zeros(size(sig_ss)); 
    rho_ss = exp(-1/chi_0); 
    q = t_0*3e-5;
    for ii = 1:6
        sig_T(ii) = sig_ss(ii)/(log(T_p/T(ii))-log(log(sqrt(rho_ss)/q))); 
        mu_T(ii) = sig_T(ii)/sqrt(rho_ss); 
    end

    for ii=7:10
        mu_T(ii) =mu_T(3);
    end
    
    r_param = mu_T(3)/G_298; 
    
    for i=1:10
        param = struct('mu_T',mu_T(i),'Kx',Kx(i),'K_pf',K_pf,'t_0',t_0,'T_p',T_p,'chi_0',0.25,'r',r_param,'nu',0.381); 
        y0=[0.0 1e-3 0.167];
        sol = ode15s(@(t,y)kclee(t,y,T(i),srate(i),param),espan{i},y0);
        ei=e{i};
        sigi=sig{i};
        ei=reshape(ei,1,[]);
        sigi=reshape(sigi,1,[]);
        err(i)=mean((deval(sol,ei,1)-sigi).^2)/(sig_ss(i))^2;
    end    
    object=sum(err);