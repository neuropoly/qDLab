function [out MSDR MHIND]=scd_model_NOGSE_full(Ax)
% [out out1 out2]=scd_model_NOGSE_full(x,G,N,R,opt,WM_param)
% Example: opt.TE=100,WM_param.fh=0.48,[out out1 out2]=scd_model_NOGSE_full(linspace(1,15,15),0.040E-3,6,4,opt,WM_param)
% x: vector containing the range of values for the CPMG part of the NOGSE
% sequence 
% G: integer, gradient value in mT/m
% N: integer, number of oscillations in the NOGSE sequence
% R: integer, axon diameter
%
% opt.verbose: switch, to plot figures or not
% opt.Color: string, indicates the color of the plotted data
% opt.onediam: switch, indicates a single diameter value or distribution
%
% WM_param.fh: hindered fraction of the intracellular volum (arbitrary units)
% WM_param.Dh: Diffusion coefficient for the hindered compartement (um^2/ms)
% WM_param.Dr: Diffusion coefficient for the restricted compartement (um^2/ms)
% WM_param.fscf:

x=Ax.x; G=Ax.G; N=Ax.N; onediam=Ax.onediam; R=Ax.WM_param.R; Dh=Ax.WM_param.Dh; Dr=Ax.WM_param.Dr; fh=Ax.WM_param.fh; var=Ax.WM_param.var;

if isfield(Ax,'verbose'), verbose=Ax.opt.verbose; else verbose =0; Ax.opt.verbose=verbose; end
if ~strcmp(Ax.norm.method,'fit'), Ax.WM_param.norm=1; end
if strcmp(Ax.norm.method,'max'), Ax.WM_param.norm=Ax.norm.maxvalue*max(Ax.data); end

y=linspace(max(x)*max(N),max(x),length(x));
if isfield(Ax,'y'); y=Ax.y; end
if ~isfield(Ax,'plotfit'), Ax.plotfit=0; end
% if Ax.Xcheck
%     xmin=(G/Ax.SR); xmin=xmin/sqrt(2);
%     xmax=Ax.opt.TE/N;
%     x= [0 (linspace(0,1,Ax.NbPts-1)).^2*(xmax-xmin)+xmin];
%     y=Ax.opt.TE-(N-1)*x;
% end

% if isfield(WM_param,'fh'); fh=WM_param.fh; else fh=0; end
% if isfield(WM_param,'Dh'); Dh=WM_param.Dh; else Dh=1.2; end
% if isfield(WM_param,'Dr'); Dr=WM_param.Dr; else Dr=0.7; end
% if isfield(WM_param,'fcsf'); fcsf=WM_param.fcsf; else fcsf=0; end

gamma=2*pi*42.57; %kHz/mT
x=x(:);y=y(:);
%stddiam=sqrt(var);
plotrand=randn*Ax.plotfit;
% gamma distribution
resol = 0.3;
%var = stddiam^2; 
b=var/R; a = R/b;
if onediam, diam = R; Er_coeff=1; else diam=0.1:resol:10;  % um : axonal diameter range
    w=gampdf(diam,a,b);
    Er_coeff = w.*(pi*diam.^2)./(pi*sum(diam.^2.*w*resol));
end

%==========RESTRICTED SIGNAL=========
MSDR=zeros(length(x),1);
for i_diam=1:length(diam)
    MSDR=MSDR+Er_coeff(i_diam)*scd_model_NOGSE(Ax);
end
if (verbose),figure (52),
    title('Restricted signal as a function x','FontSize',16)
    xlabel('x in ms (Duration of a CPMG pulse)','FontSize',14)
    ylabel('Amplitude of the normalized signal','FontSize',14)
    hold on
    plot(x,MSDR); 
end

%==========Hindered SIGNAL=========
MHIND=exp(-1/12*gamma^2*G.^2.*Dh.*((N-1)*x.^3+y.^3));
if (verbose),
    figure (6)
    hold on
    title('Hindered signal as a function x','FontSize',16)
    xlabel('x in ms (Duration of a CPMG pulse)','FontSize',14)
    ylabel('Amplitude of the normalized signal','FontSize',14)
    hold on
    plot(x,MHIND);
end

%==========FULL SIGNAL=========
Mfull=fh*MHIND+(1-fh)*MSDR;
if plotrand>1
    figure(58)
    clf
    hold on
    plot(x,Mfull,'r');
	ylim([0 1.2])
    title('Overall signal as a function of x','FontSize',16)
    xlabel('x in ms (Duration of a CPMG pulse)','FontSize',14)
    ylabel('Amplitude of the normalized signal','FontSize',14)
    if isfield(Ax,'data')
        plot(x,Ax.data/Ax.WM_param.norm,'b+');
        hold off
    end
end
hold off
%end
out=Ax.WM_param.norm*Mfull;
