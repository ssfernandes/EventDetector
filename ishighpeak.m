%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [isevent,peak]=ishighpeak(x,alpha)
%------------------------------
% INPUT
%   x [double array]: vector of values to search
%   alpha [int]: level of abnormality to search either 1.5 (standard
%       outlier) or 3 (extreme outlier (default is 3)
%------------------------------
% OUTPUT
%  isevent [bool]: true if there is an isolated high peak
%  peak [int]: index of the high peak
%------------------------------
% DESCRIPTION
%   The function searches for isolated high outliers in x
%------------------------------

if nargin==1
    alpha=3;
end

n=length(x);
[peaks, locs]=sort(x,'descend');
[~,growth_indexes]=sort(diff(x),'descend');
[~,decay_indexes]=sort(diff(x),'ascend');
ismaxpeak=or(or((decay_indexes(1)-growth_indexes(1))==1,locs(1)==1),locs(1)==n);
   
Q=quantile(x,[0.25,0.5,0.75]); Q1=Q(1); Q3=Q(3); IQ=Q3-Q1;
   
isoutlier=and(peaks(1)>Q3+alpha*IQ,peaks(1)-peaks(2)>0.05);
peak=locs(1);    
isevent=and(isoutlier,ismaxpeak);