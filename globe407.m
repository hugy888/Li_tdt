function out=globe407(e,sig)
    dict=[500,2000,1,1e-12,15000];
    gs = GlobalSearch;
    rng default % For reproducibility
    sevenmin = @(dict)obj407(dict,e,sig);
    problem = createOptimProblem('fmincon','x0',dict,...
    'objective',sevenmin,...
    'lb',[1e1,1e3,0.5,1e-13,14000],...
    'ub',[1e4,1e4,200,5e-12,16000],'nonlcon',@nlinconst);
    d = run(gs,problem);

    out=d;