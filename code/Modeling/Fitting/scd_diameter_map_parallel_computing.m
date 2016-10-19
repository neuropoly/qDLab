function result=scd_diameter_map_parallel_computing(data_raw,scheme,varargin)
% scd_diameter_map_parallel_computing(data_normed,scheme_file,...)
% OPTIONS: (default)
% output ('Results/')
% plots ('0')
% parallel ('0')
% bootstrap ('1')
% GPD : ('1')
% norm : ('fit')
% sigma_noise=1/SNR ('0.05') --> if isvector Per slice
% fitting_param ('')

% analyse diffusion data

warning('off','MATLAB:HandleGraphics:noJVM')
set(0,'DefaultFigureWindowStyle','docked')
colordef white

p = inputParser;
addRequired(p,'data_normed');
addRequired(p,'scheme');
addOptional(p,'sigma_noise','0.05');
addOptional(p,'SNR','');
addOptional(p,'selectdirection','0',@(x) ~isempty(str2num(x)));
addOptional(p,'onediameter','1',@(x) ~isempty(str2num(x)));
addOptional(p,'output','ResultsAxCaliber',@isstr);
addOptional(p,'plots','0',@(x) ~isempty(str2num(x)));
addOptional(p,'parallel','0',@(x) ~isempty(str2num(x)));
addOptional(p,'bootstrap','1',@(x) ~isempty(str2num(x)));
addOptional(p,'GPD','1',@(x) ~isempty(str2num(x)));
addOptional(p,'Dcsf','1.5',@(x) ~isempty(str2num(x)));
addOptional(p,'Zstart','1',@(x) ~isempty(str2num(x)));
addOptional(p,'remove',[],@isnumeric);
addOptional(p,'Gmax',1e5,@isnumeric); %T/m
addOptional(p,'norm','fit', @(x) any(validatestring(x,{'fit', 'max', 'none'})));
addOptional(p,'fixDh','0',@isstr);
addOptional(p,'Dr','1.4',@isstr);
%addOptional(p,'ub',[0.99    3      10         3         0.01    90   90*ones(1,10) ],@isnumeric);
addOptional(p,'optimizer','ML', @(x) any(validatestring(x,{'ML', 'mcmc'})));
addOptional(p,'NOGSE','0',@(x) ~isempty(str2num(x)));
addOptional(p,'fitting_param','',@isstr);
addOptional(p,'mask','');
addOptional(p,'fitT2','true');
addOptional(p,'Select',[]);
addOptional(p,'optimizationfun',@scd_optimization_rician_likelihood);



parse(p,data_raw,scheme,varargin{:})
in=p.Results;

% read_data
if isstr(in.data_normed)
    [data_file, ~, ext ]= sct_tool_remove_extension(in.data_normed,1); % no extension
    data_raw = load_nii_data([data_file ext]);
else data_raw=in.data_normed; end
dims=size(data_raw);

if isstr(in.scheme), scheme = scd_schemefile_read(in.scheme); else scheme=in.scheme; end

% SEQUENCE INFO
selectdirection   =   str2num(in.selectdirection);
parallel          =   str2num(in.parallel);
bootstrap              =   str2num(in.bootstrap);
Zstart            =   str2num(in.Zstart);
opt               =   in.optimizer;
fitT2             =   in.fitT2;

onediameter       =   str2num(in.onediameter);
plots             =   str2num(in.plots);
GPD               =   str2num(in.GPD);
Dcsf              =   str2num(in.Dcsf);
fixDh             =   str2num(in.fixDh);
Dr                =   str2num(in.Dr);

if ~isempty(in.SNR)
    sigma_noise=1./load_nii_data(in.SNR);
elseif strfind(in.sigma_noise,'.nii')
    sigma_noise= load_nii_data(in.sigma_noise);
elseif isstr(in.sigma_noise)
    sigma_noise       =   str2num(in.sigma_noise)*ones(dims(1:3));
else
    sigma_noise=in.sigma_noise*ones(dims(1:3));
end
sigma_noise = double(sigma_noise);
NOGSE             =   str2num(in.NOGSE);
% Axcaliber tests
norm = in.norm; % 'fit' OR 'max'
plotfit = 0;
figures = 0;




% =========================================================================
% Don't change below
% =========================================================================
output=in.output; if ~strcmp(output(end),filesep), output=[output,filesep]; end
if ~exist(output,'dir'); mkdir(output); end


% DownSample
rm =[];
remove=in.remove;
for irm=1:2:length(remove),
    rm= [rm; find(abs(scheme(:,remove(irm))-remove(irm+1))<1e-5)]; %   find(scheme_avg(:,7)~=50 | scheme_avg(:,6)==16);
end
Gmax=in.Gmax*1e-3;
rm= [rm; find(abs(scheme(:,4))>Gmax)];
scheme(rm,:)=[]; data_raw(:,:,:,rm)=[];
%==========================================================================
% Preprocessing
%==========================================================================
% Masking
if isempty(in.mask), in.mask=ones(dims(1:3)); end

% load mask
if isstr(in.mask); in.mask=load_nii_data(in.mask); in.mask = max(in.mask,[],4); end

% find X extend
Xmin=find(mean(mean(in.mask,2),3),1,'first');
Xmax=find(mean(mean(in.mask,2),3),1,'last');
Ymin=find(mean(mean(in.mask,1),3),1,'first');
Ymax=find(mean(mean(in.mask,1),3),1,'last');
%noise map
% if isempty(dir([data_file '_cor.nii*']))
%     [data] = scd_noise_estimation(scheme,data);
%     save_avw_v2(data,[data_file '_cor'],'f',[1 1 1 1],data_file,1)
%     save_avw_v2(sigma,[data_file '_sigmanoise'],'f',[1 1 1 1],data_file,1)
% else
%     data = read_avw([data_file '_cor']);
% end
% data_file = [data_file '_cor'];

% % Average same, normalize
% if normalize
%     if exist('dmri_norm.nii'), data_file='dmri_norm';
%     else
%         data_file = scd_normalize(data_file,scheme_file,param);
%     end
% end

data_avg=double(data_raw);
scheme_avg=scheme;
% [data_avg, scheme_avg, ~, ~] = scd_average_and_order_data(data_normed,scheme,param);
dims=size(data_avg);


%==========================================================================
% READ DATA
%==========================================================================

% Variables for plots
Deltavals=unique(scheme_avg(:,5));
ND=size(Deltavals,1);

if ~isempty(in.fitting_param)
    loaded_param=load(in.fitting_param); loaded_param=loaded_param.Ax(1);
    try loaded_param = rmfield(loaded_param,'scheme'); end
    try loaded_param = rmfield(loaded_param,'data'); end
    try loaded_param = rmfield(loaded_param,'plotfit'); end
    fn=fieldnames(loaded_param);
end
% Initialization of structure for parfor
for X = 1 : dims(1)
    Ax(X).data=[]; Ax(X).scheme = scheme_avg; Ax(X).fitname = '';
    for iD=1:ND
        leg{X}{iD} =  ['Delta = ' num2str(Deltavals(iD))];
        Ind_DELTA_avg{X}{iD} = find(scheme_avg(:,5)==Deltavals(iD));
    end

    Ax(X).NOGSE = NOGSE;
    Ax(X).optimizationfun = in.optimizationfun;
    if ~isempty(in.Select)
        Ax(X).Select = in.Select;
    end
    if ~isempty(in.fitting_param)
        for nf=1:length(fn), Ax(X).(fn{nf})=loaded_param.(fn{nf}); end
    else
        Ax(X).onediam = onediameter;
        Ax(X).plotfit = plotfit;
        Ax(X).figures = figures;
        Ax(X).Dcsf = Dcsf;
        Ax(X).norm = norm;
        Ax(X).fixDh     = fixDh;
        Ax(X).GPD = GPD;
        Ax(X).Dr = Dr;
        Ax(X).sigma_noise = sigma_noise;
        Ax(X).fitT2 = fitT2;
        
    end
end
Ax = orderfields(Ax);

% define q
q = scheme_avg(:,8);

% reload previous fitting



j_progress('Axonal Diameter Estimation in each voxel. Please wait...')
if parallel
    matlabpool
end



%==========================================================================
% FITTING STARTS HERE:
%==========================================================================
disp(['loop over voxels.. 0/' num2str((Xmax-Xmin+1)*(Ymax-Ymin+1)*(dims(3)-Zstart+1))])
cont=0;

% init
[X,Y,Z]=find3d(in.mask); Ax(1).data=squeeze(data_avg(X(1),Y(1),Z(1),:));
[tmp,~,pn] = Ax(1).optimizationfun(Ax(1)); Nparam=length(tmp); parametersnames = pn.parametersnames;%


% loop (in the future: reshape data in 2D (voxel,time) and reshape again at the end...)
for Z=1:dims(3)
    fitted = zeros(dims(1),dims(2),bootstrap,Nparam);
    %     if exist([output 'fitted_slice' num2str(Z) '.mat']), load([output 'fitted_slice' num2str(Z)]); else fitted = zeros(dims(1),dims(2),17); end
    
    for X = Xmin:Xmax
        % for X = 1:dims(1)
        fitted_temp = zeros(1,dims(2),bootstrap,Nparam);
        for Y= Ymin:Ymax
            if in.mask(X,Y,Z)~=0
                for ibootstrap=1:bootstrap
                    
                    AxX=Ax(X);
                    AxX.sigma_noise=AxX.sigma_noise(min(end,X),min(end,Y),min(end,Z));
                    if isfield(AxX,'Select')
                        index=find(AxX.Select);
                    else
                        index = sort(randsample(dims(4),floor(dims(4)),false));
                    end
                    %==========================================================================
                    % AxCaliber: CORE OF THE CODE !!!!!!!!!!!!
                    %==========================================================================
                    AxX.data = squeeze(data_avg(X,Y,Z,index)); AxX.scheme = AxX.scheme(index,:);
                    if strcmp(opt,'ML')
                    [fitted_temp(1,Y,ibootstrap,:), Model_signal] = AxX.optimizationfun(AxX); %
                    
                    elseif strcmp(opt,'mcmc')
                    LOOP                         = scd_fitting_MCMC(AxX.data,AxX.scheme,'plotfit',0,'perturbation',[0.01 0.01 0.1],'NN',10);
                    fitted_temp(1,Y,ibootstrap,1:17)  = [mean(LOOP.parameters) 0.5    0.01   0.5    1*ones(1,11)];                        
                    end
                    

                    % Compute SNR
                    if ~AxX.NOGSE
                        if plots
                            
                            % Plot fitting results
                            h=figure(46);
                            clf(h)
                            colordef(h,'white')
                            hold off
                            scd_display_fits(fitted_temp(1,Y,ibootstrap,:),AxX)
                            fitname = ['X' num2str(X) 'Y' num2str(Y) 'Z' num2str(Z)];
                            filename = [output, 'plots/' fitname];
                            iSaveplotX( h, filename )
                        end
                    end
                end
                
            end
            cont=cont+1;
            disp(['loop over voxels.. ' num2str(cont) '/' num2str((Xmax-Xmin+1)*(Ymax-Ymin+1)*(dims(3)-Zstart+1))])
        end
        fitted(X,:,:,:) = fitted_temp;
        
    end
    
    % save fit
    if exist([output 'fitted_results.m'])
        results=load([output 'fitted_results.m']);
    end
    results.fitted(:,:,Z,:,:)=fitted;
    clear fitted; fitted=results.fitted;
    save([output 'fitted_results'],'fitted','parametersnames')
end
disp('done')

if isfield(Ax(1),'Select')
    Select = ~~Ax(1).Select; rmfield(Ax,'Select');
    save([output 'Selected_data'],'Select')
end
save([output 'fitting_param'],'Ax')
result=scd_readresults(in.data_normed,output,onediameter, scheme_avg);


if parallel
    matlabpool close
end

end


function iSaveX( fname, x )
save( fname, 'x' );
end

function iSaveplotX( h, filename )
folder_plots=fileparts(filename);
set(gcf, 'InvertHardCopy', 'off');
if ~exist(folder_plots,'dir'), mkdir(folder_plots); end
print(h, filename,'-dtiff');
end