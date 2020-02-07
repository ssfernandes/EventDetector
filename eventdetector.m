%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [F,Tro,general_info,components,events_log]=eventdetector(T,L,F,options)
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
%   Tr [cell]: CP-APR decomposition factor matrices (Tr{i} is the factor
%       matrix associated to mode i
%   F [int]: number of relevant components
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
general_info=[];

%>> process time windows
for t=L:L:tt
    %stage 1: estimate rank of the window and decompose it
    Tro{t}{1}=[];
    Tro{t}{2}=[];
    Tro{t}{3}=[];
    
    for i=1
        [Tr,F]=decompose_window(T(:,:,t-L+1:t),F);
        Tro{t}{1}=[Tro{t}{1},Tr{1}];
        Tro{t}{2}=[Tro{t}{2},Tr{2}];
        Tro{t}{3}=[Tro{t}{3},Tr{3}];
    end

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
            if event_components(n,end)==1
                [noisy, fit, component,score]=noisecheck_highpeak(T(:,:,t-L+1:t),Tr,event_components(n,1),event_components(n,2),options);
            else
                [noisy, fit, component,score]=noisecheck_lowpeak(T(:,:,t-L+1:t),Tr,event_components(n,1),event_components(n,2),options);
            end
            
            %discard component if noisy
            if noisy
               noisy_events=[noisy_events,n];
            else
               if ~isequal(anomaly, [t,event_components(n,2)])
                   if length(identical_anomalies)==1
                       events_log=[events_log;t,F,event_components(n-1,:),prev_score];

                       ncomponents=ncomponents+1;
                       components{ncomponents}.A=identical_anomalies{1}.event{1};
                       components{ncomponents}.B=identical_anomalies{1}.event{2};
                       components{ncomponents}.C=identical_anomalies{1}.event{3};
                       components{ncomponents}.t=t;
                       
                   elseif  length(identical_anomalies)>1%select most representative anomaly in case they are high correlated
                       %discard redundant components
                       [anomalies, inds,scores,details]=findidenticanomalies(identical_anomalies);

                       for i=1:length(anomalies)
                           events_log=[events_log;t,F,details(i,:),scores(i)]; 

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
           events_log=[events_log;t,F,identical_anomalies{1}.details,prev_score];

           ncomponents=ncomponents+1;
           components{ncomponents}.A=identical_anomalies{1}.event{1};
           components{ncomponents}.B=identical_anomalies{1}.event{2};
           components{ncomponents}.C=identical_anomalies{1}.event{3};
           components{ncomponents}.t=t;
       elseif  length(identical_anomalies)>1%select most representative anomaly in case they are high correlated
           %discard redundant components
           [anomalies,inds,scores,details]=findidenticanomalies(identical_anomalies);

           for i=1:length(anomalies)
               events_log=[events_log;t,F,details(i,:),scores(i)]; 

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
    
    %stage 3: store temporary results
    save('results','t','events_log','components','Tro');
end