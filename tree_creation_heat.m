% Creation function of the tree structure without pruning for the heat
% equation
% Input:
%   A struct PDE with necessary parameters
% Outputs:
%   nodes: an nx * sum(lengths) array of all nodes of the tree, each column
%          represents one node.
%   lengths: a 1 * nt array where each column represents the number of
%            nodes at each time step.
%   adjacency_list: an na * sum(lengths(1:end-1)) array where the column 
%                     numbers refer to node numbers, and each column 
%                     represents the number of the adjacent to each node

function [nodes,lengths,adjacency_list]=tree_creation_heat(PDE)

% Initialize the parameters
x0 = PDE.ic;
control = PDE.control;
dt = PDE.dt;
nt = PDE.nt;
na = PDE.na;
a=PDE.a;
b=PDE.b;
nodes=x0;
vold=x0;
x0dt=dt*x0;
adjacency_list=[]; % adjacency matrix of the nodes
number_nodes=1;
s=1;
v=[]; % preallocation for the nodes vector
lengths=zeros(1,nt);
lengths(1)=1;
new_adiacenza=zeros(na,1); 
contnodes=1; 
nx=length(x0);
for time=2:nt % loop over discretized time domain
  for j=1:s % loop over the nodes
      for k=1:na % loop over discretized control domain
              new=tridiag(a,b,b,vold(:,j)+x0dt*control(k),nx); % new node of the tree
              v=[v new];
              contnodes=contnodes+1; % index of the node
              new_adiacenza(k)=contnodes; 
      end
     adjacency_list=[adjacency_list new_adiacenza];
  end
  vold=v;
  s=size(v,2);
  nodes=[nodes v];
  v=[];
  lengths(time)=s;
  number_nodes=number_nodes+s;
end
end