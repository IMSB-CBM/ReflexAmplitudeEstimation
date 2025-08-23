%% Main script
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
% Main script to execute to run the simulations for Schmid et al. 
%
% An excitatory stimulus is applied at random times to a motoneuron on top of a persistent
% input to simulate monosynaptic excitation (H-reflex). 
%
% can be run as parallel code
% 
clear all, close all
%% Define Scenario -------------------------------------------------------
number_of_cores = 1; % code can be parallelized
% Inputs
NumberOfMNs = 200;
num_inputs = 200; % Number of Perturbations
CycleDuration = 1000; % ms (mean inter-stimulus interval)
dt = 0.1; % ms (time step size)
%
% Definition of the excitatory stimulus 
PSC_amp = (exp(log(100).*linspace(0,1,NumberOfMNs))./100.*(0.006 - 0.01)+0.01)'; % exp distribution of amplitude in uA
PSC_tau = 1; % Time constant in ms 
single_PSC = PSC_amp./(PSC_tau/dt).*linspace(0,40/dt,40/dt+1).*exp(1-linspace(0,40/dt,40/dt+1)./(PSC_tau/dt)); % single twitch
%
% Characterize cortical input
CI_mean = 0.004; % mean cortical input in microA
CI_CoV_com = 20; % CoV of bandpass filtered noise component
CI_CoV_ind = 0.25*CI_CoV_com; % CoV of lowpass filtered noise component
%
% Control randomization
rand_seed_CI = 0;
rand_seed_perturbation = 0;
%
%% Generate Motor Neuron Pool
% Generate an object containing an empty motor neuron pool
mn_pool = MotorNeuronPool;
% Initialize the Motor Neuron Pool with a specified number of MNs and 
% default parameters
mn_pool = mn_pool.init_motor_neuron_pool(NumberOfMNs);
% Distribute the Parameters (Cm (uF), Ri (kOhm cm), Rmd, Rms (kOhm cm²), ld, ls, rd, rs (cm))
mn_pool = mn_pool.distribute_parameters(1,0.07,[14.4 6.05],...
    [1.15 0.65],[0.55 1.06], [77.5e-4 113e-4],[20.75e-4 46.25e-4],[38.75e-4 56.5e-4]);
% mn_pool = mn_pool.distribute_parameters(1,0.07,[14.4 6.05],...
   % [1.15 0.65],[0.55 1.06], [77.5e-4 113e-4],[20.75e-4 46.25e-4],[38.75e-4 56.5e-4],'exp',10);
%
% To Do: Add non default parameters here!
%
%% Run the model
parpool(number_of_cores); % starts up parallel system
% Preallocate output vectors
mn_States = zeros(number_of_cores,size(mn_pool.InitStates,1),NumberOfMNs,CycleDuration*num_inputs/number_of_cores/dt+1);
t_stim = zeros(number_of_cores,num_inputs/number_of_cores);
PSC = zeros(number_of_cores,NumberOfMNs,CycleDuration*num_inputs/number_of_cores/dt+1);
CI = zeros(number_of_cores,NumberOfMNs,CycleDuration*num_inputs/number_of_cores/dt+1);
%
% Run the simulation
parfor n = 1:number_of_cores
[mn_States(n,:,:,:),t_stim(n,:)] = psc2spiketrain(mn_pool,dt,num_inputs/number_of_cores,CycleDuration,CI_mean,CI_CoV_com,CI_CoV_ind,...
    single_PSC,rand_seed_CI+n,rand_seed_perturbation+n) 
end
%
%% Peristimulus analysis:
% Preallocate output variables:
firing_times = cell(1,number_of_cores);
spike_trains = cell(1,number_of_cores);
%
% Run peristimulus analysis and summarize results from all cores
for n=1:number_of_cores
    
    % Calculate the required time vectors
    times = linspace(0,CycleDuration*size(t_stim,2),CycleDuration*size(t_stim,2)/dt+1);
    
    % Generate spike trains from soma membrane potential
    [~,firing_times{n},~] = pot2freq(squeeze(mn_States(n,2,:,:)),dt,times(1),30);
    
    % Call the peristimulus analysis script
    LL = 300; % ms (lower limit for peristimulus analysis)
    UL = 300; % ms (upper limit for peristimulus analysis)
    bw = 1; % ms (bin size for peristimulus analysis
	[PSTH_new,PSF_new]=PeriStimAnalysis(times(squeeze(t_stim(n,:))),LL,UL,bw,firing_times{n});
    if n ==1    % Initialize variables with values of first core
       PSF = PSF_new;
       PSTH=PSTH_new;   
    else  % Save values of additional cores into the same variables  
        for mn = 1:NumberOfMNs
            PSF{mn}=[PSF{mn},PSF_new{mn}];   % Add new discharges
            PSTH{mn}(1,:) = PSTH{mn}(1,:)+PSTH_new{mn}(1,:);% sum up events per time frame
        end
    end
    
end
%
% Create the output variable for further analysis
[eval] = reflex_amplitude(PSTH,PSF,num_inputs);
%
%% Functions
% this function executes the computation of the motoneuron responses
function [mn_States,stm_Vec] = psc2spiketrain(mn_pool,dt,num_runs,CycleDuration,CI_mean,CI_CoV_com,CI_CoV_ind,...
    single_PSC,rand_seed_CI,rand_seed_perturbation) 
%
    % Vectors with time steps
    times = 0:dt:num_runs*CycleDuration; % Time steps for MN pool
    ndt = length(times);
    %
    % Create the inputs:
    
    % Cortical input
    CI = CorticalInput(CI_mean,CI_CoV_com,length(times),1000/dt,rand_seed_CI,CI_CoV_ind,mn_pool.NumberOfMUs);    
    
    % Additional input
    PSC = zeros(mn_pool.NumberOfMUs,length(times)); % Initialize PSC vector 

    % Stimulation times
    rng(rand_seed_perturbation,'twister'); % Make sure the rand alogortihm gives the same values for every run
    reg_Vec = (CycleDuration/2)/dt:CycleDuration/dt:ndt-(CycleDuration/2)/dt; % perturbation times
    rand_Vec = -100/dt+(200/dt)*rand(1,length(reg_Vec));  % Variation in perturbation times +/- 100 ms
    stm_Vec = round(reg_Vec+rand_Vec); % stimulation times
    
    PSC(:,stm_Vec) = 1; 
    for mn = 1:mn_pool.NumberOfMUs
        temp = conv(single_PSC(mn,:),PSC(mn,:)); % apply single PSC at stimulation times
        PSC(mn,:) = temp(1:ndt);
    end

    % Solve the model
    
    % Initial conditions for first time step: Initstates 
    mn_y0 = reshape(mn_pool.InitStates,6*mn_pool.NumberOfMUs,[]);
    
    options=odeset('InitialStep',dt,'RelTol',1e-5,'AbsTol',1e-5);
    
    % Preallocate mn_States
    mn_States = zeros(6,mn_pool.NumberOfMUs,ndt);
    
    % Initialize mn States
    mn_States(:,:,1) = reshape(mn_y0,6,mn_pool.NumberOfMUs,[]); 

    for time=1:length(times)-1
        % Get the MN drive for the current time step
        mu_drive = CI(:,time)+PSC(:,time);
        % Set up the integration bounds for the ODE problem
        tspan = [times(time) times(time+1)];
	% Solve model using matlab's ode23 solver
        [t,y]=ode23(@(t,y)mn_pool.evaluate_rhs_mn_pool(t,y,mu_drive'),...
                    tspan,mn_y0,options);
        mn_y0 = y(end,:)'; % New initial condition
        mn_States(:,:,time+1) = reshape(mn_y0,6,mn_pool.NumberOfMUs,[]); 
    end
    
end