%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [noisy,fit,C,score]=noisecheck_highpeak(T,Tr,f,t,options)
%------------------------------
% INPUT
%   T [sptensor]: tensor window
%   Tr [cell]: CP-APR decomposition factor matrices (Tr{i} is the factor
%       matrix associated to mode i)
%   f [int]: number of components 
%   t [int]: event instant
%   options [struct]: struct with the options of the method
%       subgraph_size [int]: minimum subgraph size to be considered as anomaly (default 0)
%       delta [double]: minimum density of the nodes in the anomalous subgraph (default 0)
%------------------------------
% OUTPUT
%   noisy [bool]: true if the candidate passed the verification test
%   fit [double]: fit of the anomalous reconstructed subgraph with respect
%       to the original one
%   C [cell]: event component (C{i} is the factor associated associated to mode i
%   score [double]: event score (rate of intervenients)
%------------------------------
% DESCRIPTION
%   The function carries out the verification on the high peak event candidate based
%   on the number of nodes participating in it, is level of deviation from
%   the remaining instants, ...
%------------------------------

if nargin<4
    alpha=0;
    delta=0;
elseif ~ismember('subgraph_size', fields(options))
    alpha=0;
    delta=options.delta; 
elseif ~ismember('delta', fields(options))
    alpha=options.subgraph_size;
    delta=0;
else
    alpha=options.subgraph_size;
    delta=options.delta; 
end

%get window length
L=size(Tr{3},1);
N=max(size(T,1:2));


%>> pattern refinement
%1. get nodes participating in the anomaly
nodes_A=find(Tr{1}(:,f));
nodes_B=find(Tr{2}(:,f));

%2. get subnetwork
Ts=reshape(double(T(nodes_A,nodes_B,t)),[length(nodes_A),length(nodes_B)]);
Ts(Ts>0)=1;

S=reshape(double(full(ktensor(1,Tr{1}(:,f),Tr{2}(:,f),Tr{3}(t,f)))),[size(T,1),size(T,2)]);
Ss=double(S(nodes_A,nodes_B));
Ss(Ss>0)=1;

%3. discard irrelevant nodes
if length(nodes_A)<=length(nodes_B)
    excess_nodes_A=find(sum(and(Ts~=Ss, Ss==1),2)>delta*length(nodes_B));
    Tr{1}(nodes_A(excess_nodes_A),f)=0;
    nodes_A(excess_nodes_A)=[];
    
    if ~isempty(nodes_A)
        %- recompute subnetwork
        Ts=reshape(double(T(nodes_A,nodes_B,t)),[length(nodes_A),length(nodes_B)]);
        Ts(Ts>0)=1;
        Ss=double(S(nodes_A,nodes_B));
        Ss(Ss>0)=1;

        excess_nodes_B=find(sum(and(Ts~=Ss, Ss==1),1)>delta*length(nodes_A));
        Tr{2}(nodes_B(excess_nodes_B),f)=0;
        nodes_B(excess_nodes_B)=[]; 
    end
else
    excess_nodes_B=find(sum(and(Ts~=Ss, Ss==1),1)>delta*length(nodes_A));
    Tr{2}(nodes_B(excess_nodes_B),f)=0;
    nodes_B(excess_nodes_B)=[];
    
    if ~isempty(nodes_B)
        %- recompute subnetwork
        Ts=reshape(double(T(nodes_A,nodes_B,t)),[length(nodes_A),length(nodes_B)]);
        Ts(Ts>0)=1;
        Ss=double(S(nodes_A,nodes_B));
        Ss(Ss>0)=1;

        excess_nodes_A=find(sum(and(Ts~=Ss, Ss==1),2)>delta*length(nodes_B));
        Tr{1}(nodes_A(excess_nodes_A),f)=0;
        nodes_A(excess_nodes_A)=[];
    end
end

Ss=double(S(nodes_A,nodes_B));
Ss(Ss>0)=1;
S(S>0)=1;

%>>compute how much of the pattern is present in each time stamp and how much
if and(~isempty(nodes_A),~isempty(nodes_B))
    for i=1:L
       %1.select subnetwork
       Ts=reshape(double(T(nodes_A,nodes_B,i)),[length(nodes_A),length(nodes_B)]);

       anomaly_presence(i)=nnz(and(Ts>0, Ss==1))/max(nnz(Ss),1);

       density(i)=nnz(Ts)/prod(size(Ts));
       mean_degree(i)=mean([sum(Ts,1)';sum(Ts,2)]);
    end
    %if the pattern corresponds to an anomaly then it is expected to:
    %- maximize the presence of the pattern in the network
    [~,t0]=ishighpeak(anomaly_presence,1.5);
    [~,t1]=ishighpeak(density,1.5);
    [~,t2]=ishighpeak(mean_degree,1.5);

    noisy=(t1~=t || t2~=t||t0~=t);
else
    noisy=true;
end
subgraph_size=length(unique([nodes_A; nodes_B]));
noisy=or(noisy, subgraph_size<=alpha);
score=subgraph_size/N;

%get output
U=double(T(:,:,t));
U(U>0)=1;
fit=nnz(U-S)/nnz(U);
C={};
C{1}=Tr{1}(:,f);
C{2}=Tr{2}(:,f);
C{3}=Tr{3}(:,f);


 