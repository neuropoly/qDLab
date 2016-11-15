function [xopt, data_model, Ax]=S0_T2_D(Ax)
% xopt=scd_optimization_rician_likelihood(Ax)
%
% INPUT:
% Ax.scheme : protocol Nx9 : Gx Gy Gz |G|(mT/um) Delta(ms) delta(ms) TE(ms) q(um-1) identifier
% Ax.data : vector of MRI data in voxel i
% Ax.sigma_noise : level of noise
%
% OUPUT:
% xopt : fitted parameters
% data_model : Synthetic MRI signal
% Ax.parametersnames : name of the fitted parameters
%

% define user options for the GUI
if ~nargin, xopt = {'T2',true}; return; end


% Define parameters to fit
Ax.parametersnames = { 'S0','T2', 'D'};
 
[S0, T2, D] = scd_assess_S0_T2_from_b0(Ax.scheme, Ax.data, 0, 1000); 
xopt = [S0 T2 D];


% Compute Smodel
b = scd_scheme_bvalue(Ax.scheme);
data_model = S0.*exp(-Ax.scheme(:,7)/T2).*exp(-b*D);
