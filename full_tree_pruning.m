% Creation function of the tree structure with pruning
% Input:
%   A struct PDE with necessary parameters
% Outputs:
%   nodes: an nx * sum(lengths) array of all nodes of the tree, each column
%          represents one node.
%   lengths: a 1 * nt array where each column represents the number of
%            nodes at each time step.
%   adjacency_list: an na * sum(lengths(1:end-1)) array where the column 
%                     numbers refer to node numbers, and each column 
%                     represents the nodes adjacent to each node


function [nodes,lengths,adjacency_list]=full_tree_pruning(PDE)

% Initialize the parameters
x0 = PDE.ic;
control = PDE.control;
dt = PDE.dt;
nt = PDE.nt;
na = PDE.na;
tol = PDE.tol;
a=PDE.a;
b=PDE.b;

% Initialize the nodes
nodes=x0;
nodes=single(nodes);
vold=x0;
x0dt=dt*x0;
nx=length(x0);
adjacency_list=[];
s=1;
v=[];
lengths=zeros(1,nt);
lengths(1)=1;
tolinv=1/tol;
Psi=single(PDE.psi*tolinv); % Divide the direction of maximum variance by the pruning tolerence
flo=floor(Psi*x0); % Assign integer indices for intervals('buckets') of length `tol` % COME BACK (maybe explaination email)
M=containers.Map(flo,[1 1]); % Map to determine the nodes contained in each bucket
tol2=single(tol*tol);
new_adiacenza=zeros(na,1);
countnodes=1;
buck=0;


for time=2:nt  % loop over discretized time domain
    
    for j=1:s % loop over the nodes of the previous level
        cont=1;
        for k=1:na % loop over discretized control domain
            
            new=single(tridiag(a,b,b,vold(:,j)+x0dt*control(k),nx)); % new node
            flo=floor(Psi*new); % project the new node onto the direction of Psi
            neighbours=[]; % determine the neighbours of the new node
            if isKey(M,flo)
                buck=1;
                neighbours=M(flo);
            end
            if isKey(M,flo-1)
                neighbours=[neighbours M(flo-1)];
            end
            if isKey(M,flo+1)
                neighbours=[neighbours M(flo+1)];
            end
            flag=0;
            if ~isempty(neighbours)
                [flag,index]=check(new,nodes(:,neighbours),tol2); % check the pruning criteria
            end
            if flag % new node is pruned         
                new_adiacenza(cont)=neighbours(index);
            else % new node is added to the tree
                v=[v new];
                countnodes=countnodes+1;
                nodes(:,countnodes)=new;
                new_adiacenza(cont)=countnodes;
                
                if buck
                    M(flo)=[M(flo) countnodes]; % update the map of buckets
                else
                    M(flo)=countnodes;
                end
            end
            buck=0;
            cont=cont+1;
        end
        adjacency_list=[adjacency_list new_adiacenza];
    end
    vold=v;
    s=size(v,2);
    v=[];
    lengths(time)=s;
end
end