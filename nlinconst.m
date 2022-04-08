function [c,ceq] = nlinconst(x)
    t_0=x(4);
    T_p=x(5);
    q = t_0*3e-5;
    c(1)=398*(-log(q)-1/(2*0.25))-T_p;
    ceq = [];
end