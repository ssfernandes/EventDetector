%---------------------------------------------------
%AUTHORS: Sofia Fernandes, Hadi Fanaee-T, Joao Gama
%--------------------------------------------------

function [event_components,nevents]=detectevents(X)
%------------------------------
% INPUT
%   X [double]: matrix of size FxT whose rows correspond to a temporal factor
%------------------------------
% OUTPUT
%   event_components [double]: matrix whose i^th row contains the details
%       of the i^th event found (namely, its index in X, the instant of time 
%       within the time window in which the event occurred and the type of
%       event (1 - highpeak; 2 - lowpeak)
%   nevents [int]: number of events detected
%------------------------------
% DESCRIPTION
%   The function searches for extreme outliers in the rows of X
%------------------------------

F=size(X,1);
event_components=[];
nevents=0;

for i=1:F
  type=1;
  [highpeak,point]=ishighpeak(X(i,:)); 

  if highpeak
      event_components=[event_components;i,point];
      nevents=nevents+1;
  end
end
