function test_ode_int()
%%% unit tests for time stepper in Vlasov 1dES Lagrangian Particle code

    delt = .5;
    tlist = 0:delt:2*pi;
    Nt = length(tlist);
    
    method_params = struct('method','rk2','delt',delt,'periodic',0,'a',1);
    ode_params = struct('function','first_order_linear_ode','lambda',1);
%     xlist = zeros([1,Nt]);
%     xlist(1) = 1;
%     for ii = 2:Nt
%         xlist(ii) = ode_int(xlist(ii-1),ode_params,method_params);
%     end
%     
%     figure
%     plot(tlist,exp(tlist),tlist,xlist)
%     legend('analytic','computed')
%     
%     sum(abs(xlist-exp(tlist)))
    
    
        ode_params = struct('function','second_order_imaginary','a',1);
    xlist = zeros([2,Nt]);
    xlist(:,1) = [1;0];
    for ii = 2:length(tlist)
        xlist(:,ii) = ode_int(xlist(:,ii-1),ode_params,method_params);
    end
    
    figure
    plot(tlist,xlist)
    figure
    plot(xlist(1,:),xlist(2,:))
    hold on
    plot(cos(tlist),-sin(tlist))
    hold off
    legend('computed','analytic')
    axis square
end


