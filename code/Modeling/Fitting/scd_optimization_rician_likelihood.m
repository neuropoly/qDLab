function [xopt, data_model, Ax]=scd_optimization_rician_likelihood(Ax)
% xopt=scd_optimization_rician_likelihood(Ax)
% MANDATORY:
% Ax.scheme : protocole
% Ax.data : données expériences
% Ax.sigma_noise

if ~isfield(Ax,'plotfit'), Ax.plotfit=0; end
if ~isfield(Ax,'fitT2'), Ax.fitT2=0; end
if ~isfield(Ax,'NOGSE'), Ax.NOGSE=0; end
Ax.output_signal=1;

[fields{1:5}]=deal('fh', 'Dh', 'diameter', 'diameter_std', 'fcsf');
[fields{17:20}]=deal('residus', 'noise_std', 'SNR_max', 'SNR_min');


Ax.parametersnames = {'fh','Dh','diameter_mean','diameter_STD','fcsf','',''};
seq_values=unique(Ax.scheme(:,9));
for iseq=1:length(seq_values)
    ind=find(Ax.scheme(:,9)==seq_values(iseq),1,'first');
    Ax.parametersnames{7+iseq}=['norm_TE' num2str(floor(Ax.scheme(ind,7))) 'D' num2str(floor(Ax.scheme(ind,5))) 'd' num2str(floor(Ax.scheme(ind,6)))];
end

if Ax.onediam, Ax.parametersnames{4}=''; end


%                              [fh     Dh        mean      std        fcsf    theta     Phi           S_b0 for each TE   ]
if ~isfield(Ax,'x0'),    Ax.x0=[[0.2   0.5         6        0.5       0.01      0.1         0.1           ones(1,length(unique(Ax.scheme(:,9))))       ]' ...
                                [0.4   0.6         6        0.5       0.01      0.1         0.1           ones(1,length(unique(Ax.scheme(:,9))))       ]' ...
                                [0.6  0.7         6        0.5       0.01      0.1         0.1           ones(1,length(unique(Ax.scheme(:,9))))       ]' ...
                                [0.8   0.9         6        0.5       0.01      0.1         0.1           ones(1,length(unique(Ax.scheme(:,9))))       ]' ...
                                [0.99   1         6        0.5       0.01      0.1         0.1           ones(1,length(unique(Ax.scheme(:,9))))       ]'];  end % starting point for the optimization (multiple points for different CSF penetration)
if ~isfield(Ax,'lb'),    Ax.lb=[0       0.3         5        0.1        0        0         0              0.7*ones(1,length(unique(Ax.scheme(:,9))))       ]';  end % lower bounds
if ~isfield(Ax,'ub'),    Ax.ub=[1       3         10         2        0.02     90        90              1.3*ones(1,length(unique(Ax.scheme(:,9))))       ]';  end % upper bound
A = zeros([size(Ax.x0,1) 1])'; A(3) = -0.95; A(4) = 1; b =-0.49;
if isfield(Ax,'Dcsf') && Ax.Dcsf, Ax.ub(5)=1; Ax.x0(5)=0.2; else Ax.parametersnames{5} = ''; end

% Find b=0 normalization
if Ax.fitT2
    Ax.norm = 'none'; 
    [Ax.S0, Ax.T2, Dh] = scd_assess_S0_T2_from_b0(Ax.scheme, Ax.data, 0, 1000); 
else
    Ax.S0 = scd_preproc_getIb0(Ax.data,Ax.scheme);
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

% find the best initialization
Ax.corrobj=1;
for istart = 1:size(Ax.x0,2), cost(istart) = objectivefunc(Ax.x0(:,istart),Ax); end; [~,I]=min(cost); Ax.x0 = Ax.x0(:,I);
%% FITTING
for i=1:2 % alternate between b=0 and microstructural fitting (otherwise bad conditionning)
    % fix b=0
    fixedparam=[false(1,5) 2+true(1,length(unique(Ax.scheme(:,9))))];
    if strcmp(Ax.norm,'norm'), Ax.corrobj=1; else Ax.corrobj=0; end
    [xopt, residue]=fmincon(@(x) objectivefunc(addfixparameters(Ax.x0,x,fixedparam),Ax), double(Ax.x0(~fixedparam)), [], [], [],[],double(Ax.lb(~fixedparam)),double(Ax.ub(~fixedparam)),[],optimoptions('fmincon','MaxIter',20,'display','off'));
    Ax.x0(~fixedparam)=xopt; xopt = Ax.x0;
    % find b=0
    Ax.corrobj=0;
    fixedparam=~[false(1,7) 1 1 false(1,length(unique(Ax.scheme(:,9)))-2)];
    [xopt, residue]=fmincon(@(x) objectivefunc(addfixparameters(Ax.x0,x,fixedparam),Ax), double(Ax.x0(~fixedparam)), A(~fixedparam), b, [],[],double(Ax.lb(~fixedparam)),double(Ax.ub(~fixedparam)),[],optimoptions('fmincon','MaxIter',20,'display','off','DiffMinChange',0.03));
    Ax.x0(~fixedparam)=xopt; xopt = Ax.x0;
end
    % fix b=0
    fixedparam=[false(1,7) true(1,length(unique(Ax.scheme(:,9))))];
    Ax.corrobj=0;
    [xopt, residue]=fmincon(@(x) objectivefunc(addfixparameters(Ax.x0,x,fixedparam),Ax), double(Ax.x0(~fixedparam)), [], [], [],[],double(Ax.lb(~fixedparam)),double(Ax.ub(~fixedparam)),[],optimoptions('fmincon','MaxIter',20,'display','off'));
    Ax.x0(~fixedparam)=xopt; xopt = Ax.x0;

% fit all
%[xopt, residue]=fmincon(@(x) objectivefunc(x,Ax), double(Ax.x0), A, b, [],[],double(Ax.lb),double(Ax.ub),[],optimoptions('fmincon','MaxIter',1,'display','off'));


%% OUTPUTS
if Ax.fitT2
    S0 = xopt(8)*abs(Ax.S0); T2 = xopt(9)*abs(Ax.T2);
    xopt(9) = T2;  xopt(8) = S0;
    Ax.parametersnames{9}='T2';  Ax.parametersnames{8} = 'S0'; [Ax.parametersnames{10:end}]=deal('');
    data_model = S0*exp(-Ax.scheme(:,7)./T2).*scd_model(xopt,Ax);
elseif strcmp(Ax.norm,'fit')
    data_model = scd_preproc_getIb0(Ax.data,Ax.scheme).*scd_model(xopt,Ax);
else
    data_model = scd_preproc_getIb0(Ax.data,Ax.scheme).*scd_model(xopt,Ax);
    [Ax.parametersnames{7:end}]=deal('');
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

function val=objectivefunc(x,Ax)
x(isnan(x))=1;
if ~Ax.NOGSE
    data_exp = max(1e-3,Ax.data);
    if Ax.fitT2
        S0 = x(8)*abs(Ax.S0); T2 = x(9)*abs(Ax.T2);
        data_model = S0*exp(-Ax.scheme(:,7)./T2).*scd_model(x,Ax);
    else
        data_model = Ax.S0.*scd_model(x,Ax);
    end
    %     figure(4)
    %     hold on
    %     l=get(gca,'xlim');
    if Ax.corrobj
        val=double(std(scd_model_likelihood_rician(data_exp,data_model, Ax.sigma_noise,0)));
        %
        %     val=double(std(data_model-data_exp));
    else
        val=double(-2*sum(scd_model_likelihood_rician(data_exp,data_model, Ax.sigma_noise,0)));
    end
    %     plot(l(2)+1,val,'+')
    if Ax.plotfit && randn>1
        figure(3), scd_display_fits(data_exp,data_model,Ax.scheme); drawnow
    end
    
else
    Ax.N=max(Ax.scheme(:,8)); Ax.G=Ax.scheme(:,4); Ax.x=Ax.scheme(:,5); Ax.y=Ax.scheme(:,6); Ax.TE=Ax.scheme(:,7);
    Ax.WM_param.fh=x(1); Ax.WM_param.Dh=x(2); Ax.WM_param.R=x(3); Ax.WM_param.var=x(4); Ax.WM_param.norm=x(7); Ax.WM_param.Dr=Ax.Dr;
    val=-2*sum(scd_model_likelihood_rician(Ax.data,scd_model_NOGSE_full(Ax), Ax.sigma_noise));
    
end
val(isnan(val))=1e10;
val(isinf(val))=1e10;

function x0 = addfixparameters(x0,x,fixedparam)
x0(~fixedparam)=x;