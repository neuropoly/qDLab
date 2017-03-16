function [xopt, data_model, Ax]=scd_optimization_rician_likelihood(Ax)
% xopt=scd_optimization_rician_likelihood(Ax)
% MANDATORY:
% Ax.scheme : protocole
% Ax.data : données expériences
% Ax.sigma_noise

if ~isfield(Ax,'plotfit'), Ax.plotfit=0; end
if ~isfield(Ax,'fitT2'), Ax.fitT2=0; end
if ~isfield(Ax,'NOGSE'), Ax.NOGSE=0; end
if ~isfield(Ax,'onediam'), Ax.onediam=1; end

Ax.output_signal=1;
Ax.data = double(max(1e-3,Ax.data));

Ax.parametersnames = {'fh','Dh','diameter_mean','diameter_STD','fcsf','',''};
TE_values=unique(Ax.scheme(:,7));
for iseq=1:length(TE_values)
    ind=find(Ax.scheme(:,7)==TE_values(iseq),1,'first');
    Ax.parametersnames{7+iseq}=['S0_TE' num2str(floor(Ax.scheme(ind,7)))];
end

if Ax.onediam, Ax.parametersnames{4}=''; end
if Ax.onediam<0, Ax.parametersnames{4}=''; Ax.parametersnames{3}=''; end



%                              [fh     Dh        mean      std        fcsf    theta     Phi           S_b0 for each TE   ]
if ~isfield(Ax,'x0'),    Ax.x0=[0.6     0.7        6        0.5       0.2      0.1         0.1           ones(1,length(unique(Ax.scheme(:,7))))       ]'; end
if ~isfield(Ax,'lb'),    Ax.lb=[0       0.3        3        0.1      0        0         0              0.7*ones(1,length(unique(Ax.scheme(:,7))))       ]';  end % lower bounds
if ~isfield(Ax,'ub'),    Ax.ub=[1       3         10         2        1     90        90              1.3*ones(1,length(unique(Ax.scheme(:,7))))       ]';  end % upper bound
A = zeros([size(Ax.x0,1) 1])'; A(3) = -0.95; A(4) = 1; b =-0.49;
if ~isfield(Ax,'Dcsf') || ~Ax.Dcsf, Ax.parametersnames{5} = ''; Ax.x0(5)=0; end


% Find b=0 normalization
if Ax.fitT2
    Ax.norm = 'none'; 
    [Ax.S0, Ax.T2] = scd_assess_S0_T2_from_b0(Ax.scheme, Ax.data, 0, 1000); 
    Ax.parametersnames{9}='T2';  Ax.parametersnames{8} = 'S0'; [Ax.parametersnames{10:end}]=deal('');
    Ax.lb(9)=0.2; Ax.ub(9)=2;
elseif strcmp(Ax.norm,'fit')
    Ax.S0 = scd_preproc_getIb0(Ax.data,Ax.scheme);
else
    Ax.S0 = scd_preproc_getIb0(Ax.data,Ax.scheme);
    [Ax.parametersnames{7:end}]=deal('');
end

if isfield(Ax,'noisepervoxel'), 
    scheme=Ax.scheme;
    scheme(scheme(:,4)==0,[5 6])=0;
    % find images that where repeated
    [~,c,ind]=consolidator(scheme(:,1:8),[],'count');
    cmax = max(c); % find images repeated more than 5 times (for relevant STD)
    if cmax<5, error('<strong>Your dataset doesn''t have 5 repeated measures (same bvec/bvals) --> you can''t estimate noise STD voxel-wise. use scd_noise_fit_histo_nii.m instead to estimate the noise STD.</strong>'); end
   
    repeated_measured = find(c==cmax);
    for irep=1:length(repeated_measured)
        STDs(irep)=std(Ax.data(ind==repeated_measured(irep)));
    end
    Ax.sigma_noise = mean(STDs);
end

% define model
fixedparam=cellfun(@isempty,Ax.parametersnames);
modelfun = @(x,scheme) max(eps,CHARMEDGPD(addfixparameters(Ax.x0,x,fixedparam),Ax.scheme,Ax)); % see: >> doc lsqcurvefit 

% % find the best initialization
% Ax.corrobj=1;
% for istart = 1:size(Ax.x0,2), cost(istart) = objectivefunc(Ax.x0(:,istart),Ax); end; [~,I]=min(cost); Ax.x0 = Ax.x0(:,I);

%% FITTING
% initiate with Gaussian noise assumption --> more stable fitting
[xopt, residue] = lsqcurvefit(modelfun,Ax.x0(~fixedparam),Ax.scheme,double(Ax.data),double(Ax.lb(~fixedparam)),double(Ax.ub(~fixedparam)),optimoptions('lsqcurvefit','MaxIter',20,'display','off'));
Ax.x0(~fixedparam)=xopt; xopt = Ax.x0;

% use Rician noise and fix fix b=0
fixedparam(6:end)=true;
modelfun = @(x,scheme) CHARMEDGPD(addfixparameters(Ax.x0,x,fixedparam),Ax.scheme,Ax); % see: >> doc lsqcurvefit 
[xopt, residue]=fmincon(@(x) double(-2*sum(scd_model_likelihood_rician(Ax.data,modelfun(x), Ax.sigma_noise))), double(Ax.x0(~fixedparam)), [], [], [],[],double(Ax.lb(~fixedparam)),double(Ax.ub(~fixedparam)),[],optimoptions('fmincon','MaxIter',20,'display','off','DiffMinChange',0.03));
Ax.x0(~fixedparam)=xopt; xopt = Ax.x0;

%% OUTPUTS
data_model=CHARMEDGPD(xopt,Ax.scheme,Ax);
if Ax.fitT2
    xopt(9)=xopt(9)*Ax.T2;
end
xopt(end+1) = residue;
Ax.parametersnames{end+1}='residue';

xopt=xopt(~cellfun(@isempty,Ax.parametersnames));
Ax.parametersnames = Ax.parametersnames(~cellfun(@isempty,Ax.parametersnames));

%% plot
% figure(3), scd_display_fits(Ax.data,data_model,Ax.scheme); drawnow

% disp(['Fh = ', num2str(xopt(1))])
% disp(['Dh = ', num2str(xopt(2))])
% disp(['mean = ', num2str(xopt(3))])
% disp(['std = ', num2str(xopt(4))])
% disp(['fcsf = ', num2str(xopt(5))])
% disp(['theta = ', num2str(xopt(6))])
% disp(['Phi = ', num2str(xopt(7))])
% disp(['Sb_0 = ', num2str(xopt(8))])
% disp(['T2 = ', num2str(xopt(9))])

% function val=objectivefunc(x,Ax)
% x(isnan(x))=1;
% if ~Ax.NOGSE
%     
%     if Ax.fitT2
%         S0 = x(8)*abs(Ax.S0); T2 = x(9)*abs(Ax.T2);
%         data_model = S0*exp(-Ax.scheme(:,7)./T2).*scd_model(x,Ax);
%     else
%         data_model = Ax.S0.*scd_model(x,Ax);
%     end
%     %     figure(4)
%     %     hold on
%     %     l=get(gca,'xlim');
%     if Ax.corrobj
%         val=double(std(scd_model_likelihood_rician(data_exp,data_model, Ax.sigma_noise,0)));
%         %
%         %     val=double(std(data_model-data_exp));
%     else
%         val=double(-2*sum(scd_model_likelihood_rician(data_exp,data_model, Ax.sigma_noise,0)));
%     end
%     %     plot(l(2)+1,val,'+')
%     if Ax.plotfit && randn>1
%         figure(3), scd_display_fits(data_exp,data_model,Ax.scheme); drawnow
%     end
%     
% else
%     Ax.N=max(Ax.scheme(:,8)); Ax.G=Ax.scheme(:,4); Ax.x=Ax.scheme(:,5); Ax.y=Ax.scheme(:,6); Ax.TE=Ax.scheme(:,7);
%     Ax.WM_param.fh=x(1); Ax.WM_param.Dh=x(2); Ax.WM_param.R=x(3); Ax.WM_param.var=x(4); Ax.WM_param.norm=x(7); Ax.WM_param.Dr=Ax.Dr;
%     val=-2*sum(scd_model_likelihood_rician(Ax.data,scd_model_NOGSE_full(Ax), Ax.sigma_noise));
%     
% end
% val(isnan(val))=1e10;
% val(isinf(val))=1e10;

function x0 = addfixparameters(x0,x,fixedparam)
x0(~fixedparam)=x;

function data_model=CHARMEDGPD(x,scheme,Ax)
% S0
if Ax.fitT2
    S0 = x(8)*abs(Ax.S0); T2 = x(9)*abs(Ax.T2);
    S0 = S0*exp(-Ax.scheme(:,7)./T2);
else
    S0 = Ax.S0;
end
% CHARMED
data_model = S0.*scd_model_CHARMED(x,Ax);

if Ax.plotfit && randn>1
    figure(3), scd_display_fits(Ax.data,data_model,Ax.scheme); drawnow
end
