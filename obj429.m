function object=obj429(dict,e,sig)
    espan = [0 25]; 
    T=298;
    srate=[4e-05 3e-04 3e-03 2e-02];
    sig_ss=[0.48 0.62 0.97 1.45];
    chi_0=0.25;
    G_298 = 2.83e3; 
    c0 = dict(1); 
    c1 = dict(2); 
    Kx = c0*exp(T/c1);
    K_pf=dict(3);
    t_0=1e-12;
    T_p=18500; 
    rho_ss = exp(-1/chi_0); 
    sig_T1 = 0.6/(log(T_p/298)-log(log(sqrt(rho_ss)/(t_0*3e-5)))); 
    mu_T1= sig_T1/sqrt(rho_ss); 

    r_param = mu_T1/G_298; 

    for i=1:4
        param = struct('mu_T',mu_T1,'Kx',Kx,'K_pf',K_pf,'t_0',t_0,'T_p',T_p,'chi_0',chi_0,'r',r_param,'nu',0.381); 
        y0=[0.0 7e-3 0.2];
        sol = ode15s(@(t,y)kclee(t,y,T,srate(i),param),espan,y0);
        ei=e{i};
        sigi=sig{i};
        ei=reshape(ei,1,[]);
        sigi=reshape(sigi,1,[]);
        err(i)=mean((deval(sol,ei,1)-sigi).^2)/(sig_ss(i))^2;
    end    
    object=sum(err);