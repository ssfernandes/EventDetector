%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [Tro,components,events_log]=eventdetector(T,L,F,options)
%------------------------------
% INPUT
%   T [sptensor]: tensor window
%   L [int]: window length 
%   F [int]: number of components to consider
%   options [struct]: struct with the options of the method
%       subgraph_size [int]: minimum subgraph size to be considered as anomaly (default 0)
%       delta [double]: minimum density of the nodes in the anomalous subgraph (default 0)
%------------------------------
% OUTPUT
%   Tro [cell]: CP-APR decomposition factor matrices (Tro{i} is the decomposition result of the i^th time window)
%   components [cell]: decomposition factors, corresponding to patterns, as follows:
%       component{i}.A - mode 1 (nodes) factor vector of event i
%       component{i}.B - mode 2 (nodes) factor vector of event i
%       component{i}.C - mode 3 (time) factor vector of event i
%       components{i}.t - instant of time (\tau) within the time window associated with the event
%   events_log [matrix]: matrix of size nX3 (with n being the no of events detected). It is organized as follows:
%       events_log[i,1] - last timestamp of the time window in which event i was detected 
%       events_log[i,2] - timestamp within the time window in which event i occurred  
%       events_log[i,3] - activity score associate with event i
%------------------------------
% DESCRIPTION
%   The function applies the proposed event detector
%------------------------------

%-------------------------------
%0. Set default parmeters
if nargin<4
   options.subgraph_size=0;
   options.delta=0.3;
elseif ~ismember('subgraph_size', fields(options))
    options.subgraph_size=0;
elseif ~ismember('delta', fields(options))
    options.delta=0.3;
end

%-------------------------------
%1.Pre-process data
size_T=size(T);

%2.Process the data in a sliding window with no overlapp
%>> initialize structures
tt=size_T(end);
components={};
ncomponents=0;
events_log=[];

%>> process time windows
window=1;
for t=L:L:tt
    %stage 1: estimate rank of the window and decompose it
    Tro{window}{1}=[];
    Tro{window}{2}=[];
    Tro{window}{3}=[];

    [Tr,F]=decompose_window(T(:,:,t-L+1:t),F);
    Tro{window}{1}=Tr{1};
    Tro{window}{2}=Tr{2};
    Tro{window}{3}=Tr{3};

    %stage 2: detect anomalous temporal evolution patterns
    [event_components,nevents]=detectevents(Tr{3}');
  
    if nevents>0
        %stage 3: filter events detected
        noisy_events=[];
        identical_anomalies={};
        ctr=1;
        anomaly=[0 0];
        for n=1:nevents
            %check if event corresponds to noise component
            [noisy, fit, component,score]=noisecheck_highpeak(T(:,:,t-L+1:t),Tr,event_components(n,1),event_components(n,2),options);
            
            %discard component if noisy
            if noisy
               noisy_events=[noisy_events,n];
            else
               if ~isequal(anomaly, [t,event_components(n,2)])
                   if length(identical_anomalies)==1
                       events_log=[events_log;t,event_components(n-1,2),prev_score];

                       ncomponents=ncomponents+1;
                       components{ncomponents}.A=identical_anomalies{1}.event{1};
                       components{ncomponents}.B=identical_anomalies{1}.event{2};
                       components{ncomponents}.C=identical_anomalies{1}.event{3};
                       components{ncomponents}.t=t;
                       
                   elseif  length(identical_anomalies)>1%select most representative anomaly in case they are high correlated
                       %discard redundant components
                       [anomalies, inds,scores,details]=findidenticanomalies(identical_anomalies);

                       for i=1:length(anomalies)
                           events_log=[events_log;t,details(i,2),scores(i)]; 

                           ncomponents=ncomponents+1;
                           components{ncomponents}.A=anomalies{i}{1};
                           components{ncomponents}.B=anomalies{i}{2};
                           components{ncomponents}.C=anomalies{i}{3};
                           components{ncomponents}.t=t;
                       end
                   end
                   %update anomaly log
                   anomaly=[t,event_components(n,2)];
                   identical_anomalies={};
                   identical_anomalies{1}.event=component;
                   identical_anomalies{1}.fit=fit;
                   identical_anomalies{1}.score=score;
                   identical_anomalies{1}.details=event_components(n,:);
                   ctr=2;
               else
                   %store anomaly
                   identical_anomalies{ctr}.event=component;
                   identical_anomalies{ctr}.fit=fit;
                   identical_anomalies{ctr}.score=score;
                   identical_anomalies{ctr}.details=event_components(n,:);
                   ctr=ctr+1;
               end
               prev_score=score;
            end
        end
              
        if length(identical_anomalies)==1
           events_log=[events_log;t,identical_anomalies{1}.details(2),prev_score];

           ncomponents=ncomponents+1;
           components{ncomponents}.A=identical_anomalies{1}.event{1};
           components{ncomponents}.B=identical_anomalies{1}.event{2};
           components{ncomponents}.C=identical_anomalies{1}.event{3};
           components{ncomponents}.t=t;
       elseif  length(identical_anomalies)>1%select most representative anomaly in case they are high correlated
           %discard redundant components
           [anomalies,inds,scores,details]=findidenticanomalies(identical_anomalies);

           for i=1:length(anomalies)
               events_log=[events_log;t,details(i,2),scores(i)]; 

               ncomponents=ncomponents+1;
               components{ncomponents}.A=anomalies{i}{1};
               components{ncomponents}.B=anomalies{i}{2};
               components{ncomponents}.C=anomalies{i}{3};
               components{ncomponents}.t=t;
           end
       end

        %update components
        event_components(noisy_events,:)=[];
        nevents=nevents-length(noisy_events);
    end
    
    window = window+1;
    
    %stage 3: store temporary results
    save('results','t','events_log','components','Tro');
end