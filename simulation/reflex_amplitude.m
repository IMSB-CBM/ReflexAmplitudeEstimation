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
%% Calculate reflex amplitude from PSF and PSTH 
% This function uses PSF and PSTH to calculate the reflex amplitude of a motor neuron pool based on the cumulative sum (CUSUM)
%
% Input:
%   * PSTH: Peristimulus time histogram (1st dimension: MN, first row:
%   count, 2nd row: time)
%   * PSF: Peristimulus frequency (1st dimension: MN, first row:
%   discharge rate, 2nd row: time)
%   * num_perturbations: number of stimuli
%
% Output:
%   reflex_eval: Structure containing all the calculated results 
%
function [reflex_eval] = reflex_amplitude(PSTH,PSF,num_perturbations)

% Add possibility to only input PSTH and PSF of a single motor neuron
if ~iscell(PSTH)
    PSTH_b = PSTH;
    clear PSTH
    PSTH{1} = PSTH_b;
end

if ~iscell(PSF)
    PSF_b = PSF;
    clear PSF
    PSF{1} = PSF_b;
end
    
num_MU = max(size(PSTH)); % number of motor neurons
LL = abs(PSTH{1}(2,1)); % Evaluation time before stimulus (ms)
UL = PSTH{1}(2,end); % Evaluation time after stimulus (ms)
bw = PSTH{1}(2,2)-PSTH{1}(2,1); % Sampling time step in ms

% Preallocate variables
mean_baseline_PSF = zeros(1,num_MU);
std_baseline_PSF = zeros(1,num_MU);
CoV_baseline_PSF = zeros(1,num_MU);
CoV_baseline_ISI = zeros(1,num_MU);
mean_baseline_PSTH = zeros(1,num_MU);
reflex_amplitude_PSTH = NaN(1,num_MU);
reflex_amplitude_PSF = NaN(1,num_MU);
reflex_latency_PSTH = NaN(1,num_MU);
reflex_latency_PSF = NaN(1,num_MU);
reflex_duration_PSTH  = NaN(1,num_MU);
reflex_duration_PSF = NaN(1,num_MU);
PSTH_array=zeros(length(PSTH{1}(1,:)),num_MU); 
CUSUM_PSTH = cell(1,num_MU);
CUSUM_PSTH_diff = cell(1,num_MU);
CUSUM_PSF = cell(1,num_MU);
sorted_PSF = cell(1,num_MU);
CUSUM_PSF_filtered = cell(1,num_MU);    
CUSUM_PSF_diff = cell(1,num_MU);   
error_box_PSTH = zeros(1,num_MU);
error_box_PSF = zeros(1,num_MU);
error_box_PSTH_diff = zeros(1,num_MU);
error_box_PSF_diff = zeros(1,num_MU);

% Cycle through all motor neurons

for mn = 1:num_MU
    % clear all temporary variables
    clear idx1 idx2 idx3 idx4 idx5 idx6 idx7 idx8 idx9 
    clear t1 t2 p1 p2

    % ----- Baseline -----

    % Mean baseline PSF
    idx1 = find(PSF{mn}(2,:)<0);
    if ~isempty(idx1)
        mean_baseline_PSF(mn) = mean(PSF{mn}(1,idx1));
        std_baseline_PSF(mn) = std(PSF{mn}(1,idx1));
        CoV_baseline_PSF(mn) = std_baseline_PSF(mn)/mean_baseline_PSF(mn)*100;
        CoV_baseline_ISI(mn) = std(1./PSF{mn}(1,idx1))/mean(1./PSF{mn}(1,idx1))*100;
    end % else: preallocated with zero

    % Mean baseline PSTH
    idx2 = find(PSTH{mn}(2,:)<0);
    if ~isempty(idx2)
        mean_baseline_PSTH(mn) = mean(PSTH{mn}(1,idx2));
    end % else: preallocated with zero

    % PSTH array: columns MNs, rows counts (needed for 3D graph)
    PSTH_array(:,mn) = PSTH{mn}(1,:)';

    % ----- Cumulative sum -----

    % PSTH
    CUSUM_PSTH{mn}(1,:) = cumsum(PSTH{mn}(1,:)-mean_baseline_PSTH(mn));
    CUSUM_PSTH{mn}(1,:) = CUSUM_PSTH{mn}(1,:)./num_perturbations; % Normalize by number of stimuli
    CUSUM_PSTH{mn}(2,:) = PSTH{mn}(2,:); % Add time steps
    CUSUM_PSTH_diff{mn} = diff(CUSUM_PSTH{mn}(1,:)); % change in CUSUM ("derivative")
    
    error_box_PSTH(mn) = max(abs(CUSUM_PSTH{mn}(1,1:LL/bw-1))); % Max from baseline (+ or -)
    
    % PSF
    [t,index] = sort(PSF{mn}(2,:)); % sort PSF according to time
    sorted_PSF{mn}(1,:) = PSF{mn}(1,index);
    sorted_PSF{mn}(2,:) = t;
    CUSUM_PSF{mn} = [cumsum(sorted_PSF{mn}(1,:)-mean_baseline_PSF(mn));sorted_PSF{mn}(2,:)];
    CUSUM_PSF{mn}(1,:) = CUSUM_PSF{mn}(1,:)./num_perturbations; % Normalize by number of stimuli 

    % Filter PSF CUSUM for reflex amplitude estimation (only use last value per full millisecond) 
    CUSUM_PSF_filtered{mn} = CUSUM_PSF{mn};
    
    i = 2;
    while i <= length(CUSUM_PSF_filtered{mn}(1,:))  
       if round(CUSUM_PSF_filtered{mn}(2,i),0) == round(CUSUM_PSF_filtered{mn}(2,i-1),0)
           CUSUM_PSF_filtered{mn}(:,i-1) =[];
       else
           i = i+1;
       end
    end

    CUSUM_PSF_diff{mn} = diff(CUSUM_PSF_filtered{mn}(1,:)); % change in CUSUM ("derivative")
    
    % Determine error boxes
    idx3 = find(CUSUM_PSF_filtered{mn}(2,:)<0);
    if ~isempty(idx3)
        error_box_PSF(mn) = max(abs(CUSUM_PSF_filtered{mn}(1,1:idx3(end))));
    end

% ----- Estimation of reflex amplitude -----

    % Turning points are determined by points where the slope is larger than during
    % baseline and where slope returns to values as large as during baseline

    % Calculate error boxes from derivatives: 
    error_box_PSTH_diff(mn) = max(abs(CUSUM_PSTH_diff{mn}(1:LL/bw-1)));
    
    idx3 = find(CUSUM_PSF_filtered{mn}(2,:)<0);
    if length(idx3) > 1 % More than 1 value in baseline
        error_box_PSF_diff(mn) = max(abs(CUSUM_PSF_diff{mn}(1:idx3(end)-1)));
    end

    % From PSTH CUSUM 
    % Reflex amplitude = difference between first event with slope
    % larger than baseline and next event with slope returning to
    % baseline
    idx4 = find(CUSUM_PSTH_diff{mn}(1,:)> error_box_PSTH_diff(mn)); % Find points with larger slope than baseline
    if ~isempty(idx4)
        for i = 1:length(idx4) % Test entries for following reflex response starting with the first
            t1 = idx4(i); % first turning point
            idx6 = find(CUSUM_PSTH_diff{mn}(1,t1:end) <= error_box_PSTH_diff(mn)); % Next slope <= error box
            if ~isempty(idx6) % Slope returns to baseline
                t2 = t1+idx6(1)-1; % last point with slope > error box
            else % slope does not change anymore
                t2 = idx4(end);
            end
            % Check if reflex exceeds CUSUM error box
            if CUSUM_PSTH{mn}(1,t2) > error_box_PSTH(mn)
                reflex_amplitude_PSTH(mn) = CUSUM_PSTH{mn}(1,t2)-CUSUM_PSTH{mn}(1,t1);
                reflex_latency_PSTH(mn) = CUSUM_PSTH{mn}(2,t1); % time of first turning point
                reflex_duration_PSTH(mn) = CUSUM_PSTH{mn}(2,t2)-CUSUM_PSTH{mn}(2,t1); % time of second turning point minus first
                break % Exit for loop
            end % Else: try next idx4 entry
        end
    end

    % From PSF CUSUM 
    % Reflex amplitude = difference between first event with slope
    % larger than baseline and next event with slope returning to
    % baseline
    idx6 = find(CUSUM_PSF_diff{mn} > error_box_PSF_diff(mn)); % Find points with larger slope than baseline
    if ~isempty(idx6) 
        for i = 1:length(idx6) % Test entries for following reflex response starting with the first
            p1 = idx6(i); % first turning point
            idx7 = find(CUSUM_PSF_diff{mn}(1,p1:end) <= error_box_PSF_diff(mn)); % Next slope <= error box
            if ~isempty(idx7) % Slope returns to baseline
                p2 = p1+idx7(1)-1; % last point with slope > error box
            else % slope does not change anymore
                p2 = idx6(end);
            end
            % Check if reflex exceeds CUSUM error box
            if CUSUM_PSF_filtered{mn}(1,p2) > error_box_PSF(mn)
                reflex_amplitude_PSF(mn) = CUSUM_PSF_filtered{mn}(1,p2)-CUSUM_PSF_filtered{mn}(1,p1);
                reflex_latency_PSF(mn) = CUSUM_PSF_filtered{mn}(2,p1); % time of first turning point
                reflex_duration_PSF(mn) = CUSUM_PSF_filtered{mn}(2,p2)-CUSUM_PSF_filtered{mn}(2,p1); % time of second turning point minus first
                break % exit for loop
            end % Else: try next idx6 entry
        end
    end

    
% Store results in the output variable
reflex_eval.PSTH = PSTH; 
reflex_eval.mean_baseline_PSTH = mean_baseline_PSTH;
reflex_eval.CUSUM_PSTH = CUSUM_PSTH;
reflex_eval.error_box_PSTH = error_box_PSTH;
reflex_eval.reflex_amplitude_PSTH = reflex_amplitude_PSTH;
reflex_eval.reflex_latency_PSTH = reflex_latency_PSTH;
reflex_eval.reflex_duration_PSTH = reflex_duration_PSTH;
reflex_eval.PSF = PSF;
reflex_eval.mean_baseline_PSF = mean_baseline_PSF;
reflex_eval.CUSUM_PSF = CUSUM_PSF;
reflex_eval.CUSUM_PSF_filtered = CUSUM_PSF_filtered;
reflex_eval.error_box_PSF = error_box_PSF;
reflex_eval.reflex_amplitude_PSF = reflex_amplitude_PSF;
reflex_eval.reflex_latency_PSF = reflex_latency_PSF;
reflex_eval.reflex_duration_PSF = reflex_duration_PSF;
reflex_eval.CoV_baseline_ISI = CoV_baseline_ISI;
reflex_eval.CoV_baseline_PSF = CoV_baseline_PSF;
reflex_eval.bounds = [-LL,UL,bw];
   
end