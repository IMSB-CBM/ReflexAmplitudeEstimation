%% Cortical Input (CI) 
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
function [CI] = CorticalInput(mean_CI,CCoV,ndt,fs,rand_seed,ICoV,NumOfMNs)
% Description:
%   Creates cortical input signals for the specified number of motoneurons
%   and with the specified signal properties
% Input:
%   * mean: mean value of common cortical input in uA
%   * CCoV: Coefficient of variance of common noise in %
%   * ndt: Number of time points
%   * fs: sampling frequency in 1/s
%   * rand_seed: seed for random algorithm
%   * ICoV: Coefficient of independent noise in %
%   * NumOfMNs: Number of motor neurons
% Output:
%   * CI: Cortical input signal for each MN in uA
        
% Bandpass filtered common drive + lowpass filtered independent noise
        
% Check dimensions of mean vector
if ~any(size(mean_CI) == 1)
	error('At least one dimension of mean must be 1.')
elseif size(mean_CI,1) ==1 && size(mean_CI,2) == 1 % only one value
            
	mean_CI = mean_CI.*ones(1,ndt); % constant mean
       
elseif size(mean_CI,1) ~= 1 && size(mean_CI,2) == 1 
	mean_CI = mean_CI'; % First dimension should be 1
end % else: correct dimension
        
% Define filters
% [transfer function coefficiants] = butterworthfilter(filter order, cut off frequency, filter type)
[A0,B0] = butter(2, 100/fs*2,'low'); % low pass
[A1,B1] = butter(2, [15/fs*2 35/fs*2],'bandpass'); % band pass
      
% Initialize random algorithm
rng(rand_seed,'twister'); 
        
% Create common drive (band pass filtered, same for all motor neurons)       
gNoise = randn(1,ndt); % white noise
commonNoise = filter(A1,B1,gNoise);
        
% Scale common noise 
drive_com = commonNoise./std(commonNoise).*(CCoV./100.*mean_CI);

% Create independent noise (low pass filtered, individual for each motor neuron)
gNoise = randn(NumOfMNs,ndt); % white noise
indNoise = filter(A0,B0,gNoise'); % low pass filter
       
% Scale independent noise
drive_ind = indNoise./std(indNoise).*(ICoV./100.*mean_CI');   

% Create drive (mean + common drive + independent noise)
CI = ones(NumOfMNs,1).*(mean_CI + drive_com) + drive_ind';

end
