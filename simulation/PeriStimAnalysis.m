
% Copyright 2025 Laura Schmid 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% “Software”), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions: 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
function [PSTH, PSF]=PeriStimAnalysis(stm_Vec,LL,UL,bw,firing_times)
% This function computes the instantaneous discharge rate of each motor neuron around
% a stimulus
% Input:
%   * stm_Vec: Vector containing the stimulation times in ms
%   * LL: Lower limit for the time frame around stimulation event in ms
%   * UL: Upper limit for the time frame around stimulation event in ms
%   * bw: Bin width in ms (sampling time step for PSTH)
%   * firing_times: in ms (cell array)
% Output:
%   * PSTH: Peristimulus time histogram (1st dimension: MN, first row:
%   count, 2nd row: time)
%   * PSF: Peristimulus frequency (1st dimension: MN, first row:
%   discharge rate, 2nd row: time)

% Time frames in ms
LL=round(LL,0); % Lower limit for time frame around stimulation event
UL=round(UL,0); % Upper limit for time frame around stimulation event 
bw=round(bw,0); % Sampling time step for event count for PSTH and PSF mean 
time_frame = (-LL:bw:UL); % Time frame for event count for PSTH and PSF mean

% Determine number of motoneurons
num_MN = max(size(firing_times));

% Make sure dimensions are correct
if size(firing_times,1) == 1
    firing_times = firing_times';
end   

% Prealloctae output variables
PSF = cell(num_MN,1);
PSTH = cell(num_MN,1);

% Perform Peri-Stim-Analysis using firing times

for mn = 1:num_MN %loop over number of motor neurons
    
    % Set all temporary variables back to zero
    PSFGr=0; PSTHGr=zeros(1,length(time_frame)); mnDR=[]; PSF_times=0;

    % Find firing times within the time frame around the stimulus times
    mnDR= 1./diff(firing_times{mn}./1000);  % Calculates discharge rates for every MN
    instant=firing_times{mn}(2:end); % Discharge times from the second to the last (no DR for the first one)

    for i = 1:length(stm_Vec) % Loop over all stimulation events  
        events=[]; % clear events
        events=find(instant>=(stm_Vec(i)-LL) & instant<=(stm_Vec(i)+UL)); % Find discharge times that are within the time range "stimulation event +UL and -DL"
        if ~isempty(events) % if events is NOT empty → there is a discharge within the time frame
            PSFGr=[PSFGr,mnDR(events)]; % Write discharge rates into the PSF vector
            PSF_times=[PSF_times,(instant(events)-stm_Vec(i))]; %time of these events according to stimulus time 
        end

        % Calculate PSTH and mean PSF
        for k = 1:length(time_frame) % loop over all time steps within the time frame
            events2 = []; % clear events2
            events2=find(instant>=(stm_Vec(i)+time_frame(k)) & instant<(stm_Vec(i)+time_frame(k)+bw)); % Find discharge times within one sampling time step bw
            if ~isempty(events2) % If times2 is Not empty -> there is a discharge within the time step
                PSTHGr(k)=PSTHGr(k)+length(events2); % The number of indentified events is added to the existing number for that time step
            end
        end
    end
 
    % Create output vectors
    PSFGr(1)=[]; % No discharge rate for first time step
    PSF_times(1)=[]; % No discharge rate for first time step
    PSF{mn}=[PSFGr;PSF_times]; % Discharge frequencies over time in ms 
    PSTH{mn}=[PSTHGr;time_frame]; % Number of events per time step over time in ms
end

end