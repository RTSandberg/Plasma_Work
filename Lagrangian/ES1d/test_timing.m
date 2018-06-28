%%% test timings

for ii = 1:10
    ode_params.function = 'odef_tracer';
    ode_params.f0vec = rand([1000,1]);
    ode_params.c1 = rand();
    ode_params.c2 = rand();
    ode_params.L = rand();
    ode_params.Ntr = 0;

    x = rand([20000,1]);

    odef_tracer(x,ode_params);
end