% Test script of the tree structure with heat equation
% y_{t} = sigma * y_{xx} + y_0(x) * u(t) for (x,t)\in [0,1] x [0,T]
% y(0,t) = y(1,t) = 0 
% y(x,0) = y_0(x)
% with cost functional: 
% J =  \int_t^T \delta * ||y(s)||^2 + gamma * |u(s)|^2 ds + \delta_2||y(t)||^2

clear
close all

% Final time
PDE.T = 1; 
% Boundaries of the control domain
umin=-1; 
umax=0;
% Step size of the spatial grid
dx=1.e-2;
% Nodes of the spatial grid
x=0:dx:1;
nx=length(x);

% Diffusion coefficient
sigma=0.15;
% Step size of the time grid
dt=dx^2/(2*sigma);
tf=1;
% Initial condition
y00=@(x) -x.^2+x;
PDE.ic=y00(x)';
PDE.x=x;

% Parameters of the cost functional
q=1;r=1;
PDE.delta=q;
PDE.delta2=q;
PDE.gamma=r; 

%% COMPUTATION OF THE MAXIMUM VARIANCE
% The first singular vector gives the direction of maximum variance
% We wil use this information to speed-up the pruning step
% We use just few time steps and few controls to capture the direction
% of maximum variance

% Nodes and grids of time discretization
PDE.dt = PDE.T/10;
PDE.t_grid = 0:PDE.dt:PDE.T;
PDE.nt=length(PDE.t_grid);
% Number of controls
PDE.na=2; 
PDE.control=linspace(umin,umax,PDE.na);
% Parameters of the implicit Euler/Centered difference method
PDE.a=1+2*sigma*PDE.dt/(dx^2);
PDE.b=-sigma*PDE.dt/(dx^2);
tic
% Nodes of full tree without pruning
[nodes,~,~]=tree_creation_heat(PDE);
% POD and direction of maximum variance for the tree structure
[Psi,U,~]=svd(nodes); 
U=diag(U);
toc
PDE.psi=Psi(:,1)';

%% REAL PARAMETERS
% Reset the mesh size of time domain
PDE.dt = 0.025;
PDE.t_grid = 0:PDE.dt:PDE.T;
PDE.nt=length(PDE.t_grid);
PDE.a=1+2*sigma*PDE.dt/(dx^2);
PDE.b=-sigma*PDE.dt/(dx^2);
% Tolerance for pruning
PDE.tol=2*PDE.dt^2;  
% Computation of the uncontrolled trajectory
PDE.na=1;
PDE.control=0;
y_unc_euler=tree_creation_heat(PDE); 

% Reset to controlled case
PDE.na=20;  
PDE.control=linspace(umin,umax,PDE.na);
%% TREE CREATION FULL

% TREE CREATION WITH PRUNING
tic
[nodes,lengths,adjacency_list]=full_tree_pruning(PDE);
TSA=toc;

% % Create map for nodes
% treenodes=containers.Map(1,nodes(:,1));
% for k=2:PDE.nt
%     new=sum(lengths(1:k));
%     treenodes(k)=nodes(:,1:new);
% end
% % The item treenodes(k) contains all nodes as an array before and on the
% % k-th fold of the tree

 %% RESOLUTION HJB
 % Approximate the value function with the TSA

 % Final data
 new = PDE.delta2*sum(nodes.^2);
 % Create a map for the value function and inserting the final condition
 % at final time t(nt)
 V = containers.Map(PDE.nt,new);
 % Create  a map for the indices
 index_control = containers.Map(PDE.nt,zeros(1,2));
 % Precomputions
 control2dtRuu = PDE.gamma*PDE.dt*PDE.control.^2';
 len=size(nodes,2);
 deltadt=PDE.delta*PDE.dt;
 newdeltadt=deltadt*new/PDE.delta2;

 % Iterate through the tree
for i=PDE.nt-1:-1:1
    len=len-lengths(i+1);
    len2=len-lengths(i)+1;
    new2=zeros(1,len);
    new3=zeros(1,len);
% Scroll through the subtree from level 1 to level i
    for j=1:len    
        % Discrete Dynamic Programming Principle
        [app, new3(j)]=min(new(adjacency_list(:,j))+control2dtRuu');
        new2(j)=app+newdeltadt(j);
    end
    index_control(i)=new3;
    V(i)=new2;
    new=new2;
end

%% COMPUTATION OPTIMAL CONTROL AND TRAJECTORY

alfaoptimal_euler=zeros(1,PDE.nt-1);
dynamics=zeros(1,PDE.nt);
dynamics(1)=1;
s=1;
for j=1:PDE.nt-1
    vett=index_control(j);
    alfaoptimal_euler(j)=PDE.control(vett(s));
    s=adjacency_list(vett(s),s);
    dynamics(j+1)=s;
end
opt_dyn_euler=nodes(:,dynamics(:));

% Computation of the uncontrolled and controlled cost
costfunctional_euler=zeros(1,PDE.nt);
costfunctional_euler(1)=V(1);
costfunctional_unc=zeros(1,PDE.nt);
costfunctional_unc(1)= deltadt*(sum(sum(y_unc_euler(:,2:end-1).^2))+(sum(y_unc_euler(:,1).^2+y_unc_euler(:,end).^2))/2)+PDE.delta2*sum(y_unc_euler(:,end).^2);
for time=2:PDE.nt-1
    app=V(time);
    costfunctional_euler(time)=app(dynamics(time));
    costfunctional_unc(time) = deltadt*(sum(sum(y_unc_euler(:,time+1:PDE.nt-1).^2))+ (sum(y_unc_euler(:,time).^2+y_unc_euler(:,end).^2))/2)+PDE.delta2*sum(y_unc_euler(:,end).^2);
end

app=V(PDE.nt);
costfunctional_euler(PDE.nt)=app(dynamics(PDE.nt));
costfunctional_unc(PDE.nt) =PDE.delta2*sum(y_unc_euler(:,end).^2);

%% PLOTS

plot(PDE.t_grid,costfunctional_unc,'r',PDE.t_grid,costfunctional_euler,'b','LineWidth',2)
grid
legend('Uncontrolled','Controlled')