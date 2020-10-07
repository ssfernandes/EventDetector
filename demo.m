%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------
clear all;
clc;

%load dataset
load('datasets\synthi')
T=sptensor(T);

%run detector on the whole network using a sliding window
dataset='synthi';
WL=[10,15,20]; %window length
Fs=[15,25,35]; %CP number of components
for F=Fs
    for L=WL
        [Tro,events_patterns,events_log]=eventdetector(T,L,F);
        save(strcat('results/',dataset,'_WL', num2str(L),'_F',num2str(F),'.mat'),'Tro','events_patterns','events_log')
    end
end

%get anomalous ranking
[top_instants]= top_events('results')

