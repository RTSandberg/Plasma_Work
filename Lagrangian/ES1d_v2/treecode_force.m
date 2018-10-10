% function force = forceeval(particles,weights,N0)

% need particles
% need weights
N=10000;
particles = [1 .2 0.4 .3 0.8]';
particles = 100*rand(N,1);

weights = [2 1 3 12 11]';
weights = 100*rand(N,1);
N0 = 10;
tic
tree = m671b_build_tree(particles,N0,[0,100]);
%%
force = zeros(size(particles));
node = 1; % set to root

for i = 1:length(particles)
    done = 0;
    node = 1;
    while ~done
        if isempty(tree(node).children) 
            force(i) = force(i) + sign(particles(i)-particles(tree(node).members))'...
                *weights(tree(node).members);
            done = 1;
        elseif particles(i) > tree(node).midpoint %&& tree(node).children(1)
                if isempty(tree(node).weights_left)
                    if tree(node).children(1)
                        tree(node).weights_left = ...
                            sum(weights(tree(tree(node).children(1)).members));
                    else
                        tree(node).weights_left = 0;
                    end
                end
            force(i) = force(i) + tree(node).weights_left;
            
            node = tree(node).children(2);
        else
%             if tree(node).children(2)
                if isempty(tree(node).weights_right)
                    if tree(node).children(2)
                        tree(node).weights_right = ...
                            sum(weights(tree(tree(node).children(2)).members));
                    else
                        tree(node).weights_right = 0;
                    end
                end
            force(i) = force(i) - tree(node).weights_right;
            node = tree(node).children(1);
%             end
        end
    end
end 
tree_time = toc
    %%
    tic
forcealt = zeros(size(particles));
for ii = 1:length(particles)
%     xi = xtvec(ii);
    forcealt(ii) = sign(particles(ii)-particles)'*weights;
end
direct_time = toc
l1err = sum(abs(force-forcealt));