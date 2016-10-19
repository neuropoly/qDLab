function [LOOP] = scd_fitting_MCMC(data,scheme,varargin)
%% Output of function is LOOP structure including :

%% Denomination of options :
colordef black
close all

p = inputParser;
addRequired(p,'data',@isnumeric);
addRequired(p,'scheme',@isstr);
addOptional(p,'NN',2000,@isnumeric);
addOptional(p,'perturbation',[0.1 0.1 0.1]);
addOptional(p,'plotfit',1,@isnumeric);

parse(p,data,scheme,varargin{:});
NN = p.Results.NN;
perturbation = p.Results.perturbation;

% Importing data and scheme file :
% scheme      = scd_schemefile_read(schemefile);
data_signal = abs(data(:)) ;

% Parameters defined by schemefile :
Ax = scd_scheme2struct(scheme) ;
Ax.data=data_signal ;

%% Initializing Markov Chain Monte Carlo Loop
j_progress('start MCMC on data...')

% Noise on the data :
Ax.sigma_noise = 0.055;

% Initializing the LOOP parameters :
param0 = [0.8 0.5 3] ;
LOOP.welcome         = 0 ;
LOOP.Signal(1,:)     = scd_model([param0   0.5    0.01   0.5    1*ones(1,6)],Ax);

% in colums : min and max -- in lines : parameters
paramlimits = [ 0 0 1 ; 0.999 3 15 ] ;

%% MCMC LOOP : Fitting MRI data
LOOP = scd_mcmc_fit(@(param) scd_compute_proba(Ax,param), param0, perturbation, paramlimits, NN);

%% Returning result
LOOP.Ax            = Ax ;
LOOP.Ax.data       = data_signal ;
LOOP.data          = data_signal ;
LOOP.perturbation  = perturbation;
LOOP.mcmc          = mean(LOOP.parameters);

end