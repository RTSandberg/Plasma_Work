function xlist = gaussian(num_rand,vthermal)
% generates a *hopefully* correlation-free thermal distribution
    width = 2.;
    xlist = zeros([num_rand,1]);
    for count = 1:num_rand
        y = 0; z = 1;
        while y < z 
            x = (2.0*rand -1.)*width;
            y = exp(-(x/vthermal)^2); % is there an over .5 in the exponent?
            z = rand;
        end
        xlist(count) = x;
    end