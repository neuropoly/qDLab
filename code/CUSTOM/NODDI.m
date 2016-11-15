function [xopt, data_model, Ax]=NODDI(Ax)
% xopt=scd_optimization_NODDI(Ax)
% THIS IS THE CLASSIC THREE STEP FITTING FOR NODDI
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
% %
% 
% dbstop if error
% 

% return options if no input
if ~nargin, xopt = {'model name',{'WatsonSHStickTortIsoV_B0','WatsonSHStickTortIsoVIsoDot_B0'}}; return; end


Ax.data = double(max(eps,Ax.data));

% GENERATE MODEL
if exist('MakeModel.m','file') ~= 2, errordlg('Please add the NODDI Toolbox to your Matlab Path: http://www.nitrc.org/projects/noddi_toolbox','NODDI is not installed properly'); return; end
Ax.model = MakeModel(Ax.modelName);

% isoIdx=GetParameterIndex(Ax.model.name,'b0');
% Ax.model.GD.fixed(isoIdx)=0; % gradient descent
% Ax.model.GS.fixed(isoIdx)=0; % grid search

isoIdx=GetParameterIndex(Ax.model.name,'diso');
Ax.model.GD.fixed(isoIdx)=1; % gradient descent
Ax.model.GS.fixed(isoIdx)=1; % grid search
Ax.model.GS.fixedvals(isoIdx)=Ax.Dcsf*1e-9;

Ax.protocol = SchemeToProtocol2(Ax.scheme);


Ax.parametersnames = Ax.model.paramsStr; 

[xopt] = ThreeStageFittingVoxel(double(Ax.data), Ax.protocol, Ax.model);

% SYNTHETIC DATA
scale = GetScalingFactors(Ax.model.name);
if (strcmp(Ax.model.name, 'ExCrossingCylSingleRadGPD') ||...
    strcmp(Ax.model.name, 'ExCrossingCylSingleRadIsoDotTortIsoV_GPD_B0'))
    xsc = xopt(1:(end-4))./scale(1:(end-1));
    theta = [xopt(end-3) xopt(end-1)]';
    phi = [xopt(end-2) xopt(end)]';
    fibredir = [cos(phi).*sin(theta) sin(phi).*sin(theta) cos(theta)]';
else
    xsc = xopt(1:(end-2))./scale(1:(end-1));
    theta = xopt(end-1);
    phi = xopt(end);
    fibredir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)]';
end
constants.roots_cyl = BesselJ_RootsCyl(30);


data_model = SynthMeas(Ax.model.name, xsc, Ax.protocol, fibredir, constants);

% OUTPUT
Ax.parametersnames{end+1} = 'ODI';
xopt(end+1) = atan2(1, xopt(3)*10)*2/pi;
