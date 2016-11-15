function [xopt, data_model, Ax]=TIME_DEPENDENT_DIFFUSION(Ax)
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

% options for the GUI
if ~nargin, xopt = {}; return; end

% Define parameters to fit
Ax.parametersnames = { 'fh','Dh', 'lc'};
x0                 = [ 0.3 , 2, 2 ];
lb                 = [ 0 , 0, 0 ];
ub                 = [ 1 , 3, 10 ];
% Define model
fun = @(x,scheme) insertyourmodelhere(x,scheme,Ax); % see: >> doc lsqcurvefit 

% Fit
xopt = lsqcurvefit(fun,x0,Ax.scheme,double(Ax.data),lb,ub);

% Compute Smodel
data_model = fun(xopt,Ax.scheme);
xopt(3) = sqrt(xopt(3))/0.2; % lc = A/0.2 : experimental


function Smodel=insertyourmodelhere(x,scheme,Ax)
Smodel=scd_preproc_getIb0(Ax.data,scheme).*TIME_DEPENDENT(x,Ax);

% figure(6)
% scd_display_qspacedata(Ax.data,Ax.scheme)
% hold on
% scd_display_qspacedata(Smodel,Ax.scheme,0,'x','-')

function Smodel = TIME_DEPENDENT(x,Ax)
% Fieremans, E., Burcaw, L.M., Lee, H.-H., Lemberskiy, G., Veraart, J., Novikov, D.S., 2016. In vivo observation and biophysical interpretation of time-dependent diffusion in human white matter. Neuroimage 129, 414?427.
% BURCAW_2015_longpulse
fext = x(1);
Dextinf = x(2); % bulk diffusivity
A = x(3); % Proportional to correlation length

t = Ax.scheme(:,5);
delta = Ax.scheme(:,6);

term1 = Dextinf;
term2 = (A./(2*delta.^2.*(t-delta/3)));
term3 = (t.^2.*log((t.^2-delta.^2)./t.^2) + delta.^2.*log((t.^2-delta.^2)./delta.^2) + 2*t.*delta.*log((t+delta)./(t-delta)));
Dh_perp = term1+term2.*term3;

Smodel = (1-fext) + fext*exp(-scd_scheme_bvalue(Ax.scheme).*Dh_perp);

