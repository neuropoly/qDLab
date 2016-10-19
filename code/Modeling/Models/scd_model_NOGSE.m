function MSDR=scd_model_NOGSE(Ax)
%

Dr=Ax.WM_param.Dr; R=Ax.WM_param.R; x=Ax.x; G=Ax.G; N=Ax.N; verbose=Ax.opt.verbose;

gamma=2*pi*42.57; %kHz/mT
y=Ax.y;
x=x(:); y=y(:);
% Dr=0.7; % um2/ms
tauc=0.26^2*R^2/Dr;
Dw_SE2=gamma^2*G.^2*Dr*tauc;

%===============HAHN==============
Mhahn=exp(-Dw_SE2.*tauc.*y.*(1-tauc./y.*(3+exp(-y/tauc)-4*exp(-y/(2*tauc)))));
if verbose, figure(55)
    plot(x,Mhahn)
    ylim([0 1])
    title('HAHN','FontSize',16)
    xlabel('x in ms (Duration of a CPMG pulse)','FontSize',14)
    ylabel('Amplitude of the normalized signal','FontSize',14)
end
%===============CPMG==============
A=(2*(N-1)+1)-(-1)^(N-1)*exp(-(N-1)*x/tauc);
%B=4*(-1)^(N-1)*(exp(-(N-1)*x/tauc).*(exp(-3/2*x/tauc)+exp(-1/2*x/tauc)+exp(-x/tauc))+exp(-3/2*x/tauc)+exp(-1/2*x/tauc)+exp(-2*x/tauc)*(N-1)+exp(-x/tauc)*(N-2))./(exp(-x/tauc)+1).^2;
B=4*((-1)^(N)*(exp(-(N-1)*x/tauc).*(exp(-3/2*x/tauc)+exp(-1/2*x/tauc)-exp(-x/tauc))+exp(-3/2*x/tauc)+exp(-1/2*x/tauc)+exp(-2*x/tauc)*(N-1)+exp(-x/tauc)*(N-2))./(exp(-x/tauc)+1).^2);
Mcpmg=min(exp(-Dw_SE2.*tauc.*((N-1)*x-tauc*(A+B))),1);
if verbose, figure(56)
    plot(x,Mcpmg)
    ylim([0 1])
    title('CPMG','FontSize',16)
    xlabel('x in ms (Duration of a CPMG pulse)','FontSize',14)
    ylabel('Amplitude of the normalized signal','FontSize',14)
end


%============CROSS TERM===========
E=1+exp(-y/tauc)-2*exp(-1/2*y/tauc)-2*exp(1/2*2*(x-y)/tauc)+exp((x-y)/tauc)+4*exp(1/2*(x-y)/tauc)-2*exp(1/2*(x-2*y)/tauc)-2*exp(1/2*x/tauc)-2*exp(1/2*x/tauc)+exp(x/tauc);
F=exp(-(x*N-x+y)/tauc)-2*exp(-1/2*(2*x*N-2*x+y)/tauc)-2*exp(-1/2*(-4*x+2*x*N+y)/tauc)+exp(-(-2*x+x*N+y)/tauc)+4*exp(-1/2*(-3*x+2*x*N+y)/tauc)-2*exp(-1/2*(-3*x+2*x*N+2*y)/tauc)+exp(-(N-1)*x/tauc)-2*exp(-1/2*x*(-3+2*N)/tauc)+exp(-x*(N-2)/tauc);
D=E+(-1)^N*F;
C=Dw_SE2.*tauc^2.*D./(exp(x/tauc)+1);
Mcross=min(exp(-C),1);
if verbose, figure(57)
    plot(x,Mcross)
    ylim([0 1])
    title('CROSS')
end

%==========OVERALL SIGNAL=========
MSDR=Mcpmg.*Mhahn.*Mcross;
if verbose, figure(58)
    hold on
    plot(x,MSDR);
    ylim([0 1])
    title('Overall signal as a function of x','FontSize',16)
    xlabel('x in ms (Duration of a CPMG pulse)','FontSize',14)
    ylabel('Amplitude of the normalized signal','FontSize',14)
end

if verbose, figure(59)
    hold on
    plot(x,MSDR/MSDR(1));
    title('OVERALL SIGNAL NORMALIZED')
end

