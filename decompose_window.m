%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [Tr,F]=decompose_window(T,F)
%------------------------------
% INPUT
%   T [sptensor]: tensor window
%   F [int]: number of components to consider
%------------------------------
% OUTPUT
%   Tr [cell]: CP-APR decomposition factor matrices (Tr{i} is the factor
%       matrix associated to mode i
%   F [int]: number of relevant components
%------------------------------
% DESCRIPTION
%   The function applies CP-APR to T using F components and then discards
%   the components that do not model a pattern
%------------------------------


rng('default');rng(0);

%decompose tensor
Tro=cp_apr(T,F,'maxiters',500);

%change format
Tr{1}=Tro.U{1}*diag(Tro.lambda);
Tr{2}=Tro.U{2};
Tr{3}=Tro.U{3};

%remove components with no temporal activity
irrelevant=[];
for i=1:F
   if sum(abs(Tr{3}(:,i)))==0
       irrelevant=[irrelevant,i];
   end
end
Tr{1}(:,irrelevant)=[];
Tr{2}(:,irrelevant)=[];
Tr{3}(:,irrelevant)=[];       
F=size(Tr{3},2);
