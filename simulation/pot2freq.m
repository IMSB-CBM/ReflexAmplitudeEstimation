%% Calculate firing times, frequency and spiketrain from membrane potential 
%
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
% Calculate a frequency from membrane potential of several motor neurons
% Input:
%   * pot: membrane potential (1st dimension: MN number, 2nd dimension:
%   time)
%   * dt: time step size in ms
%   * t_start: start time in ms (e.g. 0 or 0.1)
%   * threshold: spike threshold in units of mem_pot
% Output:
%   * freq: Calculated frequency in first dimension, time in second
%   dimension (no frequency for first discharge)
%   * firing_times: Spike times in ms
%   * spike_train: 1 for a spike, 0 else
%
function [freq,firing_times,spiketrain] = pot2freq(pot,dt,t_start,threshold)
    
    % Create time vector
    times = (t_start:dt:(size(pot,2)-1)*dt);
    
    % Initialize frequency cells
    freq = cell(size(pot,1),1);
    
    % Initialize firing times cell
    firing_times = cell(size(pot,1),1);
    
    % Preallocate spike train vector
    spiketrain = zeros(size(pot,1),size(pot,2));
    
    % Calculate spike train from membrane potential for each motor neuron
    for mn=1:size(pot,1)       
        [~,loc]=findpeaks(pot(mn,:),'MINPEAKHEIGHT',threshold,'MinPeakDistance',2/dt);
        spiketrain(mn,loc) = 1;
        firing_times{mn} = times(loc); % in ms
        ISI = diff(times(loc)); % in ms
        freq{mn}(2,:) = firing_times{mn}(2:end); % No ISI for first discharge
        freq{mn}(1,:) = 1./ISI.*1000;  % in Hz  
        loc = []; ISI = []; % Clear temporary variables
    end
    
end
