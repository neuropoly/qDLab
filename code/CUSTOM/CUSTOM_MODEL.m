function [xopt, data_model, Ax]=CUSTOM_MODEL(Ax)
% [xopt, data_model, Ax]=CUSTOM_MODEL(Ax)
%
% INPUT:
% Ax.scheme : protocol Nx9 : Gx Gy Gz |G|(mT/um) Delta(ms) delta(ms) TE(ms) q(um-1) id(integer)
% Ax.data : vector of MRI data in voxel i
% Ax.sigma_noise : level of noise
%
% OUPUT:
% xopt : fitted parameters
% data_model : Synthetic MRI signal
% Ax.parametersnames : name of the fitted parameters
%

dbstop in CUSTOM_MODEL.m at 20

errordlg('edit this file with your parameters and your model!!')

% Define parameters to fit
Ax.parametersnames = { 'fh','Dh'};
x0                 = [ 0.5 , 2 ];
lb                 = [ 0 ,   0];
ub                 = [ 1 ,   3];

% Define model
fun = @(x,scheme) insertyourmodelhere(x,scheme,Ax); % see: >> doc lsqcurvefit 

% Fit
xopt = lsqcurvefit(fun,x0,Ax.scheme,double(Ax.data),lb,ub);

% Compute Smodel
data_model = fun(xopt,Ax.scheme);


function Smodel=insertyourmodelhere(x,scheme,Ax)
fh=x(1);
Dh=x(2);
Smodel=scd_preproc_getIb0(Ax.data,scheme).*((1-fh)+fh*exp(-scd_scheme2bvecsbvals(scheme)*1e-3*Dh));

% Plot Fit
if isfield(Ax,'plotfit') && Ax.plotfit==1
    figure(6)
    scd_display_qspacedata(Ax.data,Ax.scheme)
    hold on
    scd_display_qspacedata(Smodel,Ax.scheme,0,'x','-')
end
