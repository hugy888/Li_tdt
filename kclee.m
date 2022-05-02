function out = kclee(t,y,T,eps_dot,param)

    mu_T = param.mu_T; 
    K_pf = param.K_pf; 
    K_x = param.Kx;
%     r = param.r; 
    chi_0 = param.chi_0; 
    T_p = param.T_p; 
    nu = param.nu; 
    G_el= param.G_el;  
    
    q_0 = eps_dot*param.t_0; 
%     G_el = mu_T/r;
    K_p = G_el/(K_pf*mu_T); 
    
    SIG = y(1); 
    RHO = y(2); 
    CHI = y(3); 

    Q = sqrt(RHO)*exp(-(T_p/T)*exp(-SIG/(mu_T*sqrt(RHO))));
    
    zeta = log(T_p/T)-log(log(sqrt(RHO)/q_0));
    rho_s = exp(-1/CHI); 

    dsig_de = 2*G_el*(1+nu)*(1-Q/q_0); 
    drho_de = K_p*SIG/(mu_T*zeta^2)*(Q/q_0)*(1-RHO/rho_s);
    dchi_de = K_x*(SIG/mu_T)*(Q/q_0)*(1-CHI/chi_0); 
    
    out = [dsig_de
           drho_de
           dchi_de];
