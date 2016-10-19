
function output=scd_model(x,Ax)

% INPUT
%==========================================================================
%   x           vector(>7)          fitted params includes: [fh    Dh   mean    std    fcsf  useless(test)   scale1  scale2  scale3  scale4  sca...]
%                                       fh: hindered diffusion fraction
%                                       Dh: diffusion coefficient for the hindered (extra-axonal) compartment
%                                       alpha and b define the gamma distribution of weights for different axon diameters (a)
%--------------------------------------------------------------------------
%   Ax			structure
%
%         MANDATORY
%
%
%           scheme        matrix(nb_dwi,9)      schemefile
%
%         OPTIONAL
%           data          vector(nb_dwi)      MRI Signal in one voxel
%           Dr            float               Diffusion restricted in um^2/ms;
%           plotfit       logical             plot and save fit
%           save_plot     logical             save final plot in .png
%           figures       int [0 5]           detail level
%           HinderedOnly  logical             Hindered Compartment only
%           onediam       logical             alpha (fitting param) = diameter
%           fitname       string              name of the save plot
%           norm.method   string              'fit' OR 'max' OR 'onefit' OR 'none'
%           norm.maxvalue float [0 1]         used with 'max' method
%           output_signal logical             output signal or diff? default = 1;

%==========================================================================
%KEEP ALL UNITS in um/mT/ms
if ~isfield(Ax,'data'), Ax.data=zeros(size(Ax.scheme,1),1); end
if isfield(Ax,'index'), index = Ax.index; else index=1:size(Ax.scheme,1); end
if isstr(Ax.scheme), Ax.scheme=scd_schemefile_read(Ax.scheme); end
bigdelta = Ax.scheme(index,5); bigdelta=bigdelta(:);
littledelta = Ax.scheme(index,6); littledelta=littledelta(:);
G = Ax.scheme(index,4); G=G(:);
Sdata=Ax.data(index); 
Sdata=Sdata(:);
x = real(x);
fh = x(1); Dh = x(2); mean_d = x(3); std_d = x(4); fcsf = x(5); % x(6) --> useless
x=[x(:)' ones(1,10)];
var = std_d^2; beta=var/mean_d; alpha = mean_d/beta;


if isfield(Ax,'Dr'), Dr = Ax.Dr; else Dr = 1.4; end % tortuosity model, D. Alexander 2008
if isfield(Ax,'plotfit'), plotfit = Ax.plotfit; else plotfit = 0; end
if isfield(Ax,'figures'), figures = Ax.figures; else figures = 0; end
if isfield(Ax,'HinderedOnly'), cmpEh = Ax.HinderedOnly; else cmpEh = 0; end
if isfield(Ax,'onediam'), onediam = Ax.onediam; else onediam = 1; end
if isfield(Ax,'fitname'), fitname = Ax.fitname; else fitname = 'fitplot'; end
if isfield(Ax,'norm'), norm = Ax.norm; else norm.method = 'fit'; end
if isfield(Ax,'output_signal'), output_signal = Ax.output_signal; else output_signal = 1; end
if ~isfield(Ax,'save_plot'), Ax.save_plot=0; end
if ~isfield(Ax,'Dcsf'), Ax.Dcsf=3; end
if ~isfield(Ax,'fixDh'), Dh=x(2); else if Ax.fixDh, Dh=Dr*fh/(1-fcsf); end; end


% if figures>0
%     disp(['fh    = ' num2str(x(1))])
%     disp(['Dh    = ' num2str(x(2))])
%     disp(['Dr    = ' num2str(Dr)])
%     disp(['Dcsf    = ' num2str(Ax.Dcsf)])
%     disp(['mean = ' num2str(mean)])
%     disp(['var     = ' num2str(x(4))])
%     disp(['fcsf  = ' num2str(x(5))])
%     disp(['scale  = ' num2str(x(7:end))])
% end


resol = 0.2;
if onediam, diam = mean_d; else diam=[0.1:resol:10]; end % um : axonal diameter range


gyro=42.58; %rad.KHz/mT
q=gyro.*Ax.scheme(:,6).*Ax.scheme(:,4); %um-1



%==========================================================================
%Index for plot
%==========================================================================

% seq = unique(Ax.scheme(:,9));
% Nb_seq = length(seq);
% seq_ind = cell(Nb_seq,1);
% % create index for each mixing time
% for iseq = 1: Nb_seq
%     seq_ind{iseq} = find(Ax.scheme(:,9)==seq(iseq));
% end


%==========================================================================
%Signal model for the hindered (extra-axonal) compartment
%==========================================================================
Eh=exp(-(2*pi*q).^2*Dh.*(bigdelta-littledelta/3));

% %==========================================================================
% %Signal model for the exchage between compartments
% %==========================================================================
% Eex=exp(-(2*pi*q).^2*Dex.*(bigdelta-littledelta/3));

%==========================================================================
%Signal model for CSF;
%==========================================================================

Ecsf=exp(-(2*pi*q).^2*Ax.Dcsf.*(bigdelta-littledelta/3));


%==========================================================================
% SIGNAL MODEL FOR INTRA-AXONAL
%==========================================================================
Er_sum=zeros(length(q),1);
% color = {'y', 'g', 'r', 'c', 'm', 'b', 'k'};

% weights for axon diameter distribution (gamma distribution)
 w=pdf('Gamma',diam,alpha,beta);
% figure(56)
% hold off
% plot(diam,w)
% drawnow;

Er_coeff = w.*(pi*diam.^2)./(pi*sum(diam.^2.*w*resol));
%Er_coeff = x(6:end); % diameter distribution model free
% plotrand = randn;
% if figures>0 && plotrand>0.3
%     if onediam
%         disp(['diameter :' num2str(diam)])
%     else
%     figure(100)
%     plot(diam,w)
%     hold on
%     plot(diam,Er_coeff,'r')
%     ylim([0 1])
%     hold off
%     end
% end

% Call analytical equations
for i_diam=1:length(diam)
    if (~onediam && Er_coeff(i_diam)>0.01) || i_diam==1
        R=diam(i_diam)/2;
        %b=(2*pi*q).^2.*(bigdelta-littledelta/3);
        %Er_sum= Er_sum + resol*Er_coeff(i_diam).*exp(-b.*scd_model_GPD_RDr(R*2,bigdelta,littledelta,Dr));
        Er_sum= Er_sum + resol*Er_coeff(i_diam).*scd_model_GPD(R,G,bigdelta,littledelta,Dr);
        % Er_sum= Er_sum + resol*Er_coeff(i_diam).*scd_model_smallpulse(R,q,bigdelta,Dr);
    end
end %a loop

% Calculating total response :
CHARMED = fh.*Eh + (1-fh-fcsf).*Er_sum./sum(resol*Er_coeff(Er_coeff>0.01)) + fcsf.*Ecsf;
CHARMED(isnan(CHARMED))=1;

if cmpEh, CHARMED = Eh; end

% Plotting results and data :
% if figures>1
%     figure(101)
%     for iseq = 1:Nb_seq
%         plot(q(seq_ind{iseq}),Er_sum(seq_ind{iseq}), 'Color',color{iseq})
%         ylim([0 1])
%         hold on
%         plot(q(seq_ind{iseq}),Eh(seq_ind{iseq}), 'b')
%     end
%     hold off
%     for iseq = 1:Nb_seq
%         figure(102)
%         plot(q(seq_ind{iseq}),Eh(seq_ind{iseq}), 'Color',color{iseq})
%         ylim([0 1])
%         hold on
%     end
%     hold off
%     
%     for iseq = 1:Nb_seq
%         figure(103)
%         plot(q(seq_ind{iseq}),Ecsf(seq_ind{iseq}), 'Color',color{iseq})
%         ylim([0 1])
%         hold on
%     end
%     hold off
%     
% end



% size_sdata=size(Sdata);
%n_Sdata=abs(Sdata)/max(abs(Sdata));
Nb_seq = length(unique(Ax.scheme(:,9)));
seqnumbering = unique(Ax.scheme(:,9));
for iseq = 1:Nb_seq
    seq_ind = Ax.scheme(:,9) == seqnumbering(iseq);
    if strcmp(norm,'fit')
        CHARMED(seq_ind)=abs(CHARMED(seq_ind))*x(iseq+7); %OR max(abs(Sdata(TM{i_Delta})))*0.85;
    elseif strcmp(norm,'max')
        Sdata(seq_ind)=abs(Sdata(seq_ind))/max(abs(Sdata))*norm.maxvalue;
        CHARMED(seq_ind)=abs(CHARMED(seq_ind))/max(abs(CHARMED(seq_ind)))*norm.maxvalue;
    elseif strcmp(norm,'onefit')
         Sdata(seq_ind)=abs(Sdata(seq_ind))/x(8);
         x(iseq+8)=x(7);
    end
end

% if plotfit
%     figure(100)
%     scd_display_qspacedata(Ax.data,Ax.scheme)
%     hold on
%     scd_display_qspacedata(CHARMED,Ax.scheme,'none','-')
%     hold off
%     drawnow;
% end


if sum(w)*resol < 0.8 && ~onediam, output = 100*ones(1,length(CHARMED));
else
    output=(CHARMED-Sdata); %./((q+0.01)/max(q));
end
if figures>0, disp(['diff = ' num2str(sum(output))]); end

if output_signal, output = CHARMED; end


