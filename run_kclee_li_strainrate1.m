clear all; %close all; 
load('exp_data.mat')

t_0 = 5e-11; 
chi_0 = 0.25;

% elastic properties
nu = 0.381; 
G_298 = 2.83e+3; 

%% STEADY STATE 
% calculation at two strain rates, T = const
T = 298; 
% q_1 = t_0*4.e-5; 
% q_2 = t_0*3e-4;
% 
% sig_1 = 0.48; 
% sig_2 = 0.62; 
rho_ss = exp(-1/chi_0); 

% A = log(log(sqrt(rho_ss)/q_2)) - log(log(sqrt(rho_ss)/q_1)); 

% sig_T = (sig_1 - sig_2)/A; 
% T_p = T*exp(sig_1/sig_T)*log(sqrt(rho_ss)/q_1);
% calculation at multiple T
q = t_0*3e-5;
% T_p=1.05*398*(-log(q)-1/(2*chi_0));
T_p=14500;
T_all = [248 273 298 348 398]; 
sig_ss = [0.99 0.71 0.6 0.49 0.37];
sig_T = zeros(size(sig_ss)); 
mu_T = zeros(size(sig_ss)); 
for ii = 1:length(T_all)
    sig_T(ii) = sig_ss(ii)/(log(T_p/T_all(ii))-log(log(sqrt(rho_ss)/q))); 
    mu_T(ii) = sig_T(ii)/sqrt(rho_ss); 
end 

% ratio between Taylor stress and shear modulus 
% r_param = mu_T(3)/G_298; 


clear A ii q q_1 q_2 sig_1 sig_2 sig_ss T 

%% calculation for T=298 and edot=3e-5
srate=[4e-05 3e-04 3e-03 2e-02];

eps_dot = 3.e-5; 
T = 298;
color=['b' 'r' 'm' 'g' 'c'];

hold on
espan = [0 0.25]; 
e = linspace(0,0.25,1000);

c0 = 1; 
c1 = 500; 
% Kx = c0*exp(T/c1);
Kx = 75;
K_pf = 1;
% K_p = 1e5; 
rho_ini = 1e-3; 
chi_ini = 0.2; 
y0 = [0.0 rho_ini chi_ini];
j=1;
legend_list={'4e-5','3e-4','3e-3','2e-2'};
for ii=1:4
    e_exp = Li(~isnan(Li(:,2*ii-1+12)),2*ii-1+12);
    s_exp = Li(~isnan(Li(:,2*ii+12)),2*ii+12);

    % mu_T(3) corresponds to T=298 K
    param = struct('mu_T',mu_T(3),'Kx',Kx,'K_pf',K_pf,'t_0',t_0,'T_p',T_p,'chi_0',chi_0,'G_el',2.83e+3,'nu',nu); 

    sol = ode15s(@(t,y)kclee(t,y,298,srate(ii),param),espan,y0);
    s = deval(sol,e,j);
    % error(4)=mean((deval(sol,e4,1)-sig4).^2);

    scatter(e_exp,s_exp,color(ii));
    plot(e*100,s,color(ii));
% xlabel('e(%)') 
% ylabel('Stress, MPa') 
end

title('Different Strain Rate')
