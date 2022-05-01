function out=globe421(e,sig)
    dict=[1,100,5];
    gs = GlobalSearch;
    rng default % For reproducibility
%     opts=optimoptions('fmincon', 'MaxFunctionEvaluations',2);
    sevenmin = @(dict)obj429(dict,e,sig);
    problem = createOptimProblem('fmincon','x0',dict,...
    'objective',sevenmin,...
    'lb',[1e-1,1,1],...
    'ub',[1e3,1e3,10]);
    d = run(gs,problem);

    out=d;