%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [components,inds,scores,details]=findidenticanomalies(A)
%------------------------------
% INPUT
%   A [cell] cell whose i^th entry is a struct containing the nodes factor
%       vector and the score for event i
%------------------------------
% OUTPUT
%   components [cell]: i^th entry consists of thecomponents (factor vectors) 
%       describing the unique patterns
%   inds [int array]:  i^th entry is the index of uniique event i in A
%   scores [double array]: i^th contains the scores of the i^th unique
%       pattern found
%------------------------------
% DESCRIPTION
%   This function compares the subgraphs associated with anomalies in A and
%   discard the anomalies which correspond to a similar event
%------------------------------


%compute correlation between events
for i=1:length(A)
    for j=i+1:length(A)
        C(i,j)=max(mean([corr(A{i}.event{1},A{j}.event{1}),corr(A{i}.event{2},A{j}.event{2})]),mean([corr(A{i}.event{1},A{j}.event{2}),corr(A{i}.event{2},A{j}.event{1})]));      
        
        if and(abs(sum(A{i}.event{1}>0)-sum(A{j}.event{1}>0))<0.3*max(sum(A{i}.event{1}>0),sum(A{j}.event{1}>0)),abs(sum(A{i}.event{2}>0)-sum(A{j}.event{2}>0))<0.3*max(sum(A{i}.event{2}>0),sum(A{j}.event{2}>0)))
            overlap(i,j)=(sum(and(A{i}.event{1}>0,A{j}.event{1}>0))+sum(and(A{i}.event{2}>0,A{j}.event{2}>0)))/(sum(A{i}.event{1}>0)+sum(A{i}.event{2}>0));
        elseif and(abs(sum(A{i}.event{1}>0)-sum(A{j}.event{2}>0))<0.3*max(sum(A{i}.event{1}>0),sum(A{j}.event{2}>0)),abs(sum(A{i}.event{2}>0)-sum(A{j}.event{1}>0))<0.3*max(sum(A{i}.event{2}>0),sum(A{j}.event{1}>0)))
            overlap(i,j)=(sum(and(A{i}.event{1}>0,A{j}.event{2}>0))+sum(and(A{i}.event{2}>0,A{j}.event{1}>0)))/(sum(A{i}.event{1}>0)+sum(A{i}.event{2}>0));
        else
            overlap(i,j)=0;
        end
        
    end
end
C(C>=0.7)=1;
C(C<1)=0;

components={};
ctr=1;
details=[];
considered=[];
redundant=[];
scores=[];
for i=size(C,2):-1:1
    if and(~ismember(i,redundant),~ismember(i,considered))
        similar=unique([find(C(:,i));find(overlap(:,i)>0.7)]);
        if isempty(similar)
            components{ctr}=A{i}.event;
            considered=[considered,i];
            scores(ctr)=A{i}.score;
            details(ctr,:)=A{i}.details;
            ctr=ctr+1;
        else
            similar=[i;similar];
            %get fit
            minfit=Inf;
            mincomp=1;
            for j=1:length(similar)
                if minfit>A{similar(j)}.fit
                    minfit=A{similar(j)}.fit;
                    mincomp=similar(j);
                end
            end

            %add component with more variance explained
            if ~ismember(mincomp, considered)
                components{ctr}=A{mincomp}.event;
                considered=[considered,mincomp];
                scores(ctr)=A{mincomp}.score;
                details(ctr,:)=A{mincomp}.details;
                ctr=ctr+1;
            end
            similar(similar==mincomp)=[];
            redundant=[redundant;similar];
        end
    end
end

inds=considered;


