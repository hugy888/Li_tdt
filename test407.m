Li = readmatrix('Li.csv');
e={};
sig={};
for i=1:10
    e_i=Li(:,2*i-1);
    e_i=e_i(~isnan(e_i));
    e=[e e_i];
    sig_i=Li(:,2*i);
    sig_i=sig_i(~isnan(sig_i));
    sig=[sig sig_i];
end
d=globe407(e,sig);
disp(d)