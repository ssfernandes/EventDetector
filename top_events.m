%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [top_instants]= top_events(path)
%------------------------------
%INPUT
% 	path [string]: location of the ensemble results
% 	T [int]: number of instants
%------------------------------
%OUTPUT
%   top_instants [array]: list of top abnormalous events  (top_instants(i) is the instant with the i^th highest abnormalous score)
%------------------------------
% DESCRIPTION
%   Tool to facilitate the reading of the ensemble results 
%------------------------------
% NOTE
% the files must be in a format <*>_WL<window_length>_<*>

%get list of results
results_files= dir(strcat(path,'/*.mat'));
results_files = cellfun(@(x) fullfile(x), {results_files.name}, 'UniformOutput', false);

%get number of timestamps in network
T=0;
for i=1:length(results_files)
    file=results_files{i};
    strcat(path,'/',file)
    load(strcat(path,'/',file));
    T=max(T,length(Tro));
end

%for each timestamp count the number of occurrences 
events_freq=zeros(T,2);
for i=1:length(results_files)
    file=results_files{i};
    load(strcat(path,'/',file));
    
    %get window length from path
    file_details=strsplit(file,'WL');
    setting=strsplit(file_details{end},'_');
    L=str2num(setting{1});
    
    %get unique events detected at the current setting
    current_events=[];
    for j=1:size(events_log,1)

        %get events
        if ~ismember(events_log(j,1)-L+events_log(j,4), current_events) && events_log(j,6)>0
            events_freq(events_log(j,1)-L+events_log(j,4),1)=events_freq(events_log(j,1)-L+events_log(j,4),1)+1;

            %store score
            if events_freq(events_log(j,1)-L+events_log(j,4),1)==1
               events_freq(events_log(j,1)-L+events_log(j,4),2)=events_log(j,6);
            end
            current_events=[current_events,events_log(j,1)-L+events_log(j,4)];       
        end
    end
end

timestamps=1:T;
[freq,rank]=sort(events_freq(:,1),'descend');
scores=events_freq(rank,2);
timestamps=timestamps(rank);

%if the number of events detected is large, filter the ones associated with
%smaller subgraphs
if sum(scores>0)>20
    th=quantile(scores(scores>0),0.5);
    freq(scores<th)=[];
    timestamps(scores<th)=[];
    scores(scores<th)=[];
else
    freq(scores==0)=[];
    timestamps(scores==0)=[];
    scores(scores==0)=[];
end

%read 1st freq
ctr=1;
tied_instants=timestamps(ctr);
tied_inds=ctr;
tied_scores=scores(ctr);
prev_freq=freq(ctr);

while ctr<length(freq)
    %move forward in the ranking
    ctr=ctr+1;
    
    %get number of models
    current_freq=freq(ctr);

    %update tied log
    if current_freq==prev_freq
        tied_instants=[tied_instants,timestamps(ctr)];
        tied_inds=[tied_inds,ctr];
        tied_scores=[tied_scores,scores(ctr)];
    else
        if length(tied_instants)>1
            [~,inds]=sort(tied_scores,'descend');
            timestamps(tied_inds)=tied_instants(inds);
        end
        tied_instants=timestamps(ctr);
        tied_inds=ctr;
        tied_scores=scores(ctr);
    end
    prev_freq=freq(ctr);
end

top_instants=timestamps;




        