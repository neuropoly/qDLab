function [xopt, data_model, Ax]=DTI(Ax)
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

Ax.data = max(double(Ax.data),eps);
% Define parameters to fit
Ax.parametersnames = { 'FA','L1','L2','L3'};

% fit
D=scd_model_dti(Ax.data./scd_preproc_getIb0(Ax.data,Ax.scheme),Ax.scheme);
[~,L]=eig(D); L = sort(diag(L),'descend');

L_mean = sum(L)/3;

% compute metrics
FA = sqrt(3/2)*sqrt(sum((L-L_mean).^2))/sqrt(sum(L.^2));

xopt = [FA, L(1), L(2), L(3)];

% Compute Smodel
bvec=Ax.scheme(:,[1 2 3]);
data_model = scd_preproc_getIb0(Ax.data,Ax.scheme).*exp(-scd_scheme_bvalue(Ax.scheme).*diag(bvec*D*bvec'));


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
