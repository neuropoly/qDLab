function [xopt, data_model, Ax]=scd_optimization_custom_model(Ax)
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

dbstop if error

error('edit this file with your parameters and your model!!')

% Define parameters to fit
Ax.parametersnames = { 'fh','Dh','diameter_mean'};
x0                 = [ 0.3 , 2  ,      4        ];
 
% Define model
fun = @(x,scheme) insertyourmodelhere(x,scheme,Ax); % see: >> doc lsqcurvefit 

% Fit
xopt = lsqcurvefit(fun,x0,Ax.scheme,double(Ax.data));

% Compute Smodel
data_model = fun(xopt,Ax.scheme);


function Smodel=insertyourmodelhere(x,scheme,Ax)
x(4:10)=[1 0 0 0 1 1 1];
Smodel=scd_preproc_getIb0(Ax.data,scheme).*scd_model(x,Ax);