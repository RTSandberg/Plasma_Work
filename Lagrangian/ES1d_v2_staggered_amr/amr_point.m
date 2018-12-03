%%% amr
%%% takes in arrays:
%%% alpha, x, v, f
%%% checks x, v for insertion condition
%%% when condition is met, inserts point
%%% returns new alpha, x, v, f, and new length
%%% 
%%% insertion condition: currently tests whether two points are far apart
%%% in phase space
%%% interpolation: linear

function [alpha, x, v, f,Ntemp] =  amr_point(alpha, x, v, f,L,delx,N,eps)
%%% additional items used: L, delx, N, 

 %%%%%%%%%amr:point insertion: test against flow, to see if need to add
    %%%%%%%%%points
    %%% use l1 norm, check to see if normalized distance in x,v space
    %%% is more than 2 cell lengths

    xsep = x(2:end) - x(1:end-1);
    xsep = min(abs(xsep),min(abs(xsep-L),abs(xsep+L)));
    vsep = abs(v(2:end) - v(1:end-1));
    testind = sqrt(xsep.^2 + vsep.^2) > eps;
    %%%%%%%%%%%%%%%%%%%%%%
    interp_method = 'cubic';
%     
    xtemp = zeros([ceil(3*N/2),1]); xtemp(1:2) = x(1:2);
    alphatemp = zeros(size(xtemp)); alphatemp(1:2) = alpha(1:2);
    vtemp = zeros(size(xtemp)); vtemp(1:2) = v(1:2);
    ftemp = zeros(size(xtemp)); ftemp(1:2) = f(1:2);
    
    Ntemp = N;
    amrcount = 2;
    
    for jj = 2:N-2
        if testind(jj)
            Ntemp = Ntemp + 1;
            alphastar = .5*(alphatemp(amrcount) + alpha(jj+1));
            alphatemp(amrcount+1) = alphastar;
            alphatemp(amrcount+2) = alpha(jj+1);
            
            if strcmp(interp_method, 'linear')
                
            fjj = ftemp(amrcount)/(alpha(jj+1)-alphatemp(amrcount-1));
            fjjp1 = f(jj+1)/(alpha(jj+2)-alphatemp(amrcount));
            
                ftemp(amrcount) = fjj*(alphastar - alphatemp(amrcount-1));
                ftemp(amrcount+1) = (fjj+fjjp1)*.5*(alpha(jj+1)-alphatemp(amrcount));
                % the conservative route is the average, linearly interpolating
                % route
                ftemp(amrcount+2) = fjjp1*(alpha(jj+2)-alphastar);

                %linear interpolation in v
                vtemp(amrcount+1:amrcount+2) = [.5*(v(jj)+v(jj+1));v(jj+1)];

                %just regular linear interpolation
                if abs(x(jj+1)-x(jj)) < L/2
                    xtemp(amrcount+1:amrcount+2) = [.5*(x(jj)+x(jj+1));x(jj+1)];
                else
                    xtemp(amrcount+1:amrcount+2) = [.5*(L+x(jj)+x(jj+1));x(jj+1)];
                end
                
            elseif strcmp(interp_method, 'cubic')
                alpha_list = [alpha(jj-1),alpha(jj),alpha(jj+1),alpha(jj+2)];
                jjm2 = mod(jj-2-1,N)+1;
                jjp3 = mod(jj+3-1,N)+1;
                fjjm1 = f(jj-1)/(alpha(jj)-alpha(jjm2));
                fjj = f(jj)/(alpha(jj+1)-alpha(jj-1));
                fjjp1 = f(jj+1)/(alpha(jj+2)-alpha(jj));
                fjjp2 = f(jj+2)/(alpha(jjp3)-alpha(jj+1));
                
                ftemp(amrcount) = fjj*(alphastar - alphatemp(amrcount-1));
                fstar = cubic_interp(alpha_list,...
                    [fjjm1,fjj,fjjp1,fjjp2]);
                ftemp(amrcount+1) = fstar*(alpha(jj+1)-alphatemp(amrcount));
                ftemp(amrcount+2) = fjjp1*(alpha(jj+2)-alphastar);
                
                vtemp(amrcount+1) = cubic_interp(alpha_list, ...
                    v(jj-1:jj+2));
                vtemp(amrcount+2) = v(jj+1);
                
                %unwrap the x
                xlist = [x(jj-1) x(jj) x(jj+1) x(jj+2)];
                for uu = 2:4
                    if xlist(uu) - xlist(uu-1) > L/2
                        xlist(uu) = xlist(uu)-L;
                    elseif xlist(uu-1) - xlist(uu) > L/2
                        xlist(uu) = xlist(uu)+L;
                    end
                end
                xtemp(amrcount+1) = mod(cubic_interp(alpha_list,xlist),L);
                xtemp(amrcount+2) = x(jj+1);
                
            end
            % an attempt at momentum conserving interpolation
            %vtemp(amrcount+1:amrcount+2) = ...
            %    [(vnew(jj)*fjj+v(jj+1)*fjjp1)/(fjj+fjjp1);v(jj+1)];
            
                
           
            amrcount = amrcount + 2;
            
        else
            alphatemp(amrcount+1) = alpha(jj+1);
            xtemp(amrcount+1) = x(jj+1);
            vtemp(amrcount+1) = v(jj+1);
            ftemp(amrcount+1) = f(jj+1);
            amrcount = amrcount + 1;
        end
%     
    end
    % since we were lazy and aren't dealing with the wraparound at the 
% endpoints, we have to manually fill in the last few points of the temp
% vectors
xtemp(amrcount+1) = x(end);
alphatemp(amrcount+1) = alpha(end);
vtemp(amrcount+1) = v(end);
ftemp(amrcount+1) = f(end);


    alpha = alphatemp(1:Ntemp);
    f = ftemp(1:Ntemp);
    v = vtemp(1:Ntemp);
    x = xtemp(1:Ntemp);
end
    
    
function ystar = cubic_interp(x,y)
xstar = .5*(x(2)+x(3));
dd12 = (y(1)-y(2)) / (x(1)-x(2));
dd23 = (y(2)-y(3)) / (x(2)-x(3));
dd34 = (y(3)-y(4)) / (x(3)-x(4));
dd123 = (dd12 - dd23) / (x(1)-x(3));
dd234 = (dd23 - dd34) / (x(2)-x(4));
dd1234 = (dd123-dd234) / (x(1)-x(4));
ystar = y(1) + (dd12 + (dd123 + dd1234*(xstar - x(3))) * (xstar - x(2)) )*(xstar - x(1));
end
   