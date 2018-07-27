function x = ode_int(x, ode_params, method_params)
%%% takes one iteration of an ode integrator
%%% input:
%%%     x:          vector, variable of ode to advance
%%%     ode_params: struct, parameters of the ode function
%%%         relevant fields:
%%%             function: string, indicates name of file with right-hand side of ODE
%%%     method_params:struct, parameters of method
%%%         all:     
%%%                 delt:   scalar, time step
%%%                 periodic:   boolean, indicates whether domain is periodic
%%%                  xmin:       scalar, needed for making method periodic
%%%                  L:          scalar, needed for making method periodic
%%%                  method:     string, options are
%%%             'fe':   forward Euler, 1st order, single step, explicit
%%%             'be':   %NOT SUPPORTED% backward Euler, 1st order, single step, implicit
%%%             'l':    %NOT SUPPORTED% staggered leapfrog, 2nd order, two step, explicit
%%%             'rk2':   Runge-Kutta, 2nd order, 
%%%                     single step, two stage, explicit
%%%             'rk4':   %NOT SUPPORTED% Runge-Kutta, 4th order, 
%%%                     single step, four stage, explicit
%%%         specific parameters:
%%%             for rk2:     a: scalar, determines which rk2 method to use
%%% output:
%%%     x:  vector, new x after one time step              

ode_function = ode_params.function;
     
method = method_params.method;
delt = method_params.delt;
periodic = method_params.periodic;
if periodic
    xmin = method_params.xmin;
    L = method_params.period;
end

if strcmp(method,'fe')
%     'Forward Euler'
    fn = eval(strcat(ode_function, '(x,ode_params)'));    
    x = x + delt*fn;
    
elseif strcmp(method,'be')
%     'implicit Euler'
elseif strcmp(method,'l')
%     'leapfrog'
elseif strcmp(method,'rk2')
%     '2nd order Runge-Kutta'
    a = method_params.a;
    c = 1/2/a; b = 1-c;

    fn = eval(strcat(ode_function, '(x,ode_params)'));
    x2 = x + a*delt*fn;
    if periodic
        N = length(x2)/2;
        x2(1:N) = mod(x2(1:N)-xmin,L)+xmin;
    end
    f2 = eval(strcat(ode_function, '(x2,ode_params)'));
    v = (b*fn + c*f2);
    
    x = x + delt*v;
elseif strcmp(method,'rk4')
%     '4th order Runge-Kutta'
    f1 = eval(strcat(ode_function, '(x,ode_params)'));
    x2 = x + delt/2*f1;
    if periodic
        N = length(x2)/2;
        x2(1:N) = mod(x2(1:N)-xmin,L)+xmin;
    end
    f2 = eval(strcat(ode_function, '(x2,ode_params)'));
    x3 = x + delt/2*f2;
    if periodic
        N = length(x3)/2;
        x3(1:N) = mod(x3(1:N)-xmin,L)+xmin;
    end
    f3 = eval(strcat(ode_function, '(x3,ode_params)'));
    x4 = x + delt*f3;
    if periodic
        N = length(x4)/2;
        x4(1:N) = mod(x4(1:N)-xmin,L)+xmin;
    end
    f4 = eval(strcat(ode_function, '(x4,ode_params)'));
    
    v = (f1 + 2*f2 + 2*f3 + f4);
    
    x = x + delt/6*v;
end

if periodic
    N = length(x)/2;
    x(1:N) = mod(x(1:N)-xmin,L)+xmin;
end