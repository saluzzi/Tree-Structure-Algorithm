% Function to deduce if the pruning rule is active

% The pruning rule is stated as following:
%    Two given nodes \zeta_i and \zeta_j can be merged if:
%    ||\zeta_i - \zeta_j || \leq \epsilon_T
%    for a given threshold \epsilon_T.

% Inputs:
%   u: an nx * 1 array representing the new node
%   nodes: an array representing the previous nodes already computed at the
%          same level
%   tol2: a given positive threshold (\epsilon_T)

% Outputs:
%   flag: flag variable such that  if flag = 1, the pruning rule is active 
%         and the new node will not be added to the tree,
%         while if flag = 0, the new node will be added to the tree.
%   index: the index of the previous node that is merged with the new
%           node. If index = 0, the new node will be added to the tree.

function [flag,index]=check(u,nodes,tol2) % flag zero if the new node is not close to any of the previous nodes
    v=sum((nodes-u).^2);
   if sum(v<tol2)>0
       [~,index]=min(v);
       flag=1;
   else
       flag=0;
       index=0;
   end
end