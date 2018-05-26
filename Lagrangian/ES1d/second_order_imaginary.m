function xprime = second_order_imaginary(x,ode_params)
    a = ode_params.a;
    xprime = [0,a; -a,0]*x;
end