function tree_out = m671b_build_tree(passed_particles,passed_N0,entire_interval) % Barnes-Hut

global tree particles node_count N0
particles = passed_particles; N0 = passed_N0;
% N = 100; N0 = 1; rand N = rand(N,1); particles = rand N(:,1);
N = length(particles); %N0 = 1; particles = [0.4186 0.8462 0.5252 0.2026 0.6721];
tree = struct( 'interval' , [] , 'members' , [] , 'children' , [] ,'midpoint',[],'weights_left',[],'weights_right',[]);
tree(1).interval = entire_interval;
tree(1).members = 1:N;
tree(1).midpoint = (entire_interval(1)+entire_interval(2))/2;
node_count = 1;
root = 1; build_tree(root,0); tree_out = tree;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function build_tree(cluster_index,single_levels)
global tree particles node_count N0
child = struct( 'interval' , [] , 'members' , [] , 'children' , [],'midpoint',[],'weights_left',[],'weights_right',[] );
n = length(tree(cluster_index).members);
if (n > N0)
%
% step 1 : define intervals for child clusters
%
 a = tree(cluster_index).interval(1); b = tree(cluster_index).interval(2);
 tree(cluster_index).midpoint = (a+b)/2;
 midpoint = tree(cluster_index).midpoint;
 child(1).interval = [a midpoint]; child(2).interval = [midpoint b];
 %
 % step 2 : insert particles from parent into child clusters
 %
 count(1) = 0; count(2) = 0;
for j = 1:n
    particle_index = tree(cluster_index).members(j);
    index = 1; 
    if particles(particle_index) > midpoint;   index = 2;    end
    child(index).members = [child(index).members particle_index];
    count(index) = count(index) + 1;
end
%
% step 3 : add non-empty children to tree
%
tree(cluster_index).children = [0 0];
for j = 1:2
  if (count(j) >= 1)
    node_count = node_count + 1;
    tree(cluster_index).children(j) =  node_count; 
    tree = [tree child(j)];
  end
end
%
% step 4 : recursive call to build next level of children
%
          
for i = 1:length(tree(cluster_index).children)
    if ~tree(cluster_index).children(1) || ~tree(cluster_index).children(2)
        single_levels = single_levels+1;
    else
        single_levels = 0;
    end
    if tree(cluster_index).children(i) && single_levels < 3
    cluster_index_new = tree(cluster_index).children(i); 
    build_tree(cluster_index_new,single_levels);
    end
end
end