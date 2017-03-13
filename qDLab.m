function varargout = qDLab(varargin)
% qDLab(dmri_fname,schemefile)
%
% qDLab MATLAB code for qDLab.fig
%      qDLab, by itself, creates a new qDLab or raises the existing
%      singleton*.
%
%      H = qDLab returns the handle to a new qDLab or the handle to
%      the existing singleton*.
%
%      qDLab('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in qDLab.M with the given input arguments.
%
%      qDLab('Property','Value',...) creates a new qDLab or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before qDLab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to qDLab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qDLab

% Last Modified by GUIDE v2.5 15-Nov-2016 12:21:18

% Begin initialization code - DO NOT EDIT
warning('off','all'); 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @qDLab_OpeningFcn, ...
    'gui_OutputFcn',  @qDLab_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before qDLab is made visible.
function qDLab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qDLab (see VARARGIN)
plotedit off
set(gcf,'WindowStyle','normal')
% Choose default command line output for qDLab
qDLabDir = fileparts(which(mfilename()));
addpath(genpath(qDLabDir));

handles.output = hObject;
handles.Z=1;
handles.X=1;
handles.Y=1;
handles.h=[];
handles.mask_fname_all = '';
if length(varargin)<1, varargin{1}=''; end
if length(varargin)<2, varargin{2}=''; end
handles = ChooseData(handles,varargin{1},varargin{2});
if ~isfield(handles,'scheme'), help scd_schemefile_create; disp('<strong>Please select a 4D NIFTI file and a schemefile.</strong>');  close(handles.figure1); return; end
handles.selecttoolindex = false([size(handles.scheme,1) 1]);

% Simplify options depending on the dataset:
    % fit T2?
if length(unique(handles.scheme(scd_scheme2bvecsbvals(handles.scheme)<1000,7))) == 1 % Don't propose to fit T2 if only one echo time
    set(handles.norm_fitT2,'Visible','off');
else
    set(handles.norm_fitT2,'Value',true);
end

    % fit diameter distribution?
%if ~sum(scd_scheme2bvecsbvals(handles.scheme)>30000), set(handles.gammadiam,'enable','off'); end % don't propose gamma distribution if maximal bvalue<30,000


    % Load Models
qDLab_dir = fileparts(mfilename('fullpath'));
CUSTOM=[qDLab_dir filesep 'code' filesep 'CUSTOM'];
addpath(CUSTOM)
CUSTOM_list = sct_tools_ls([CUSTOM filesep '*.m']);
set(handles.modelname,'String',cat(1,get(handles.modelname,'String'),CUSTOM_list'))
model_list=get(handles.modelname,'String');
if handles.qspace3D && ~isempty(find(strcmp(model_list,'NODDI.m'), 1))
    set(handles.modelname,'Value',find(strcmp(model_list,'NODDI.m'), 1))
    handles = modelname_Callback(handles.modelname,[],handles);
end
    % estimate noise voxel-wise?
[~,c]=consolidator(handles.scheme(:,1:8),[],'count');
cmax = max(c); % find images repeated more than 5 times (for relevant STD)
if cmax<5 && sum(handles.scheme(:,4)==0)<5
    set(handles.noisepervoxel,'enable','off')
elseif cmax>10 % noise using STD is robust --> select this option by default
    set(handles.noisepervoxel,'Value',true)
end

 % one slice
 if size(handles.data,3) ==1
     set(handles.Zslider,'enable','off')
 end

axes(handles.MRI)
cla;
imshow(mat2gray(mean(handles.data(:,:,handles.Z,:),4))); colormap gray;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = qDLab_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = get(handles.AcquisitionList,'UserData');


% --- Executes on selection change in AcquisitionList.
function AcquisitionList_Callback(hObject, eventdata, handles)
% hObject    handle to AcquisitionList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AcquisitionList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AcquisitionList

% --- Executes during object creation, after setting all properties.
function AcquisitionList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AcquisitionList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pixel.
function pixel_Callback(hObject, eventdata, handles)
% hObject    handle to pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pixel
axes(handles.MRI)
cla;
imshow(mat2gray(mean(handles.data(:,:,handles.Z,:),4))); colormap gray;
[handles.X,handles.Y,~]=impixel;
% color the pixel in red for the user
handles.data_pixelXY=mat2gray(repmat(squeeze(mean(handles.data(:,:,handles.Z,:),4)),[1 1 3]));
for i_point=1:length(handles.Y)
    hold on
    plot(handles.X(i_point),handles.Y(i_point),'rx');
end

set(handles.uipanel3,'Visible','on')
set(handles.uipanel1,'Visible','on')
set(handles.uipanel2,'Visible','on')
set(handles.text29,'Visible','on')





guidata(hObject, handles);
plotbutton_Callback(hObject, eventdata, handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on button press in ChooseScheme.
function handles=ChooseScheme(handles,scheme_fname_all)
% hObject    handle to ChooseScheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ChooseScheme
if isempty(scheme_fname_all)
    [scheme_fname,path] = uigetfile({'*.txt;*.scheme'},'Select schemefile');
    scheme_fname_all=[path,scheme_fname];
else
    scheme_fname = scheme_fname_all;
end
set(handles.textscheme,'String',['SchemeFile : ' scheme_fname])
handles.scheme_fname_all=scheme_fname_all;

if scheme_fname
    [handles.scheme, handles.qspace3D]=scd_schemefile_read(scheme_fname_all);
    if handles.qspace3D, set(handles.Xaxis_Gz,'Value',1);  end
    if size(handles.data,4)~=size(handles.scheme,1), error(['<strong>Error: your dataset has ' num2str(size(handles.data,4)) ' while your schemefile has ' num2str(size(handles.scheme,1)) ' rows.</strong>']); end

    % sort data
    [handles.q,Isort]=sort(handles.scheme(:,8));
    handles.data=handles.data(:,:,:,Isort); handles.scheme=handles.scheme(Isort,:);
    % list acquisitions
    [~,ia]=unique(handles.scheme(:,9),'rows');
    handles.acqList=handles.scheme(ia,[9 7:-1:5]);
        % add Gmax
    for i=1:size(handles.acqList,1), 
        Gmaxaq = max(handles.scheme(handles.scheme(:,9)==handles.scheme(ia(i),9),4))*1e6;
        handles.acqList(i,5)=Gmaxaq;
        acqListString{i}=sprintf('%5i %10.0f %10.0f %10.0f %10.0f',handles.acqList(i,:)); 
    end
    set(handles.AcquisitionList,'String',acqListString);
    set(handles.AcquisitionList,'Value',1:size(handles.acqList,1)); % select all by default
    % list directions
    handles.bvecs=unique(handles.scheme(:,1:3),'rows');
    for i=1:size(handles.bvecs,1), dirList{i}=sprintf('% 10.2f % 10.2f % 10.2f',handles.bvecs(i,:)); end
    set(handles.DirList,'String',dirList);
    set(handles.DirList,'Value',1:size(handles.bvecs,1)); % select all by default
end


% --- Executes on button press in ChooseData.
function handles = ChooseData(handles,data_fname_all,schemefile)
% hObject    handle to ChooseData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ChooseData
if isempty(data_fname_all)
    [data_fname,path] = uigetfile({'*.nii;*.nii.gz'},'Select data');handles.data_fname_all=[path,data_fname];
else
    handles.data_fname_all = data_fname_all;
    data_fname = data_fname_all;
end
if handles.data_fname_all
    set(handles.textdata,'String',['data (NIFTI) : ' data_fname])
    
    dat = load_untouch_nii(handles.data_fname_all); handles.data = dat.img;
    handles=ChooseScheme(handles,schemefile);
    % set slider extent
    if size(handles.data,3)>1
        set(handles.Zslider,'Max',size(handles.data,3))
        set(handles.Zslider,'SliderStep',[1 1]/(size(handles.data,3)-1))
        set(handles.Zslider,'Value',round(size(handles.data,3)/2))
        handles = Zslider_Callback(handles.Zslider,[],handles);
    else
        set(handles.Zslider,'HandleVisibility','off')
    end
    
    % set Gmax
    Gmax = max(handles.scheme(:,4))*1e6;
    set(handles.Gmax_edit,'String',num2str(Gmax));
    
    handles.data=handles.data;
end


% --- Executes during object creation, after setting all properties.
function textdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in plotbutton.
function plotbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotbutton
axes(handles.plot1)
legend('off')
if ~get(handles.hold,'Value')
    cla;
end

% Get Selected Acquisitions
acqselected=get(handles.AcquisitionList,'Value');
dirselected=get(handles.DirList,'Value'); bvecs=handles.bvecs(dirselected,:);

    % Directions
dir_logical=false;
for i_dir=1:length(dirselected)
    dir_logical=dir_logical | ismember(handles.scheme(:,1:3),bvecs(i_dir,:),'rows');
end

% Gmax
dir_logical = dir_logical & handles.scheme(:,4)<=(str2double(get(handles.Gmax_edit,'String'))*1e-6+eps);

    % Acquisitions
handles.Selection=false;
for i=1:length(acqselected)
    handles.Selection=handles.Selection + i.*(dir_logical & handles.scheme(:,9)==handles.acqList(acqselected(i),1) & ~handles.selecttoolindex);
end

if get(handles.jet,'Value')
    handles.colorplot=jet(length(acqselected)*length(handles.Y));
elseif get(handles.cool,'Value')
    handles.colorplot=cool(length(acqselected)*length(handles.Y));
end


S0=0;
% plot
for i_point=1:length(handles.Y) % loop on data point
    % Compute fitted models
    mdel=[];
    if isfield(handles,'modelfit') && ~isempty(handles.modelfit)
        mdel=handles.modelfit{i_point};
    end
    
    for iaq=setdiff(unique(handles.Selection),0)' % loop on acquisitions parameters
        % change color for different point
        plotnumber=iaq+max(handles.Selection)*(i_point-1);
        handles.data_pixelXY(handles.Y(i_point),handles.X(i_point),:)=handles.colorplot(plotnumber,:);
        axes(handles.MRI)
        cla;
        imagesc(handles.data_pixelXY);
        % == plot data == 
        axes(handles.plot1)
        schemeiaq = handles.scheme(handles.Selection==iaq,:);
        data = squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,handles.Selection==iaq));
        
%         Merge same qvalues:
%         [~,ia,iq]=unique(schemeiaq(:,8)); schemeiaq = schemeiaq(ia,:);
%         [~,data]=consolidator(iq,data);
        
            % abscissa
        switch get(get(handles.Abscissa,'SelectedObject'),'Tag')
            case 'Xaxis_bvalue'
                bvals=scd_scheme2bvecsbvals(schemeiaq); 
                absc=bvals; 
            case 'Xaxis_qvalue'
                qvalues=schemeiaq(:,8);
                absc=qvalues; 
            case 'Xaxis_Gz'
                [~,~,fiberdirection] = scd_model_dti(squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,:))./scd_preproc_getIb0(squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,:)),handles.scheme),handles.scheme);
                Gnorm = 1;
                Gz=schemeiaq(:,1:3)*fiberdirection(:);
                absc = Gz./Gnorm;
            case 'XaxisCustom'
                absc_all=scd_display_Xaxiscustom(handles.scheme,squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,:)));
                absc = absc_all(handles.Selection==iaq);
        end
        [~,order]=sort(absc); absc = absc(order); data = data(order); 
        if ~isempty(mdel)
            modeliaq = mdel(handles.Selection==iaq); modeliaq = modeliaq(order);
        end
            % normalize data
            if get(handles.plotting_normalize,'value')
                normvalue=scd_preproc_getIb0(squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,:)),handles.scheme);
                normvalue=normvalue(handles.Selection==iaq);
                data = data./normvalue;
            else
                normvalue = 1;
            end
            handles.h(plotnumber)=plot( absc,data,'+','Color',handles.colorplot(plotnumber,:));
            set(handles.h(plotnumber),'MarkerSize',15,'LineWidth',2)
            hold on
            
        % Don't show data button
        if ~get(handles.DisplayData,'Value')
            set(handles.h(plotnumber),'Visible','off');
        end
        % Show legend for the first point only
        if i_point==1
            set(handles.h(plotnumber),'DisplayName',['G_{max}=' num2str(handles.acqList(acqselected(iaq),5),'%.0f') 'mT/m \Delta=' num2str(handles.acqList(acqselected(iaq),4),'%.0f') 'ms \delta=' num2str(handles.acqList(acqselected(iaq),3),'%.0f') 'ms TE=' num2str(handles.acqList(acqselected(iaq),2),'%.0f') 'ms']);
        else
            hAnnotation = get(handles.h(plotnumber),'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
        end
        
        
        % plot model
        if ~isempty(mdel)

                handles.g(iaq)=plot(absc,modeliaq(:)./normvalue,'Color',handles.colorplot(plotnumber,:));
                set(handles.g(iaq),'Linewidth',3)       
                hAnnotation = get(handles.g(iaq),'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')

        end


    end
end
ylim('auto')
limits=get(gca,'Ylim'); set(gca,'Ylim',[0 limits(2)]);
xlim('auto')
guidata(hObject, handles);
% display legend
legid=legend(handles.plot1,'show');
set(legid,'FontSize',20)
set(handles.plot1,'XColor',[1 1 1])
set(handles.plot1,'YColor',[1 1 1])



% --- Executes on button press in DisplayData.
function DisplayData_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DisplayData


% --- Executes on slider movement.
function handles = Zslider_Callback(hObject, eventdata, handles)
% hObject    handle to Zslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.Zslider,'Value',floor(get(handles.Zslider,'Value')));
handles.Z=get(handles.Zslider,'Value');
guidata(hObject, handles);
if ~isempty(handles.data)
    axes(handles.MRI)
    cla;
    imshow(mat2gray(squeeze(mean(handles.data(:,:,handles.Z,:),4))))
end

% --- Executes during object creation, after setting all properties.
function Zslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in DirList.
function DirList_Callback(hObject, eventdata, handles)
% hObject    handle to DirList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DirList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DirList


% --- Executes during object creation, after setting all properties.
function DirList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hold.
function hold_Callback(hObject, eventdata, handles)
% hObject    handle to hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold


% --- Executes on button press in jet.
function jet_Callback(hObject, eventdata, handles)
% hObject    handle to jet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jet
set(handles.cool,'Value',~get(handles.cool,'Value'))

% --- Executes on button press in cool.
function cool_Callback(hObject, eventdata, handles)
% hObject    handle to cool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cool
set(handles.jet,'Value',~get(handles.jet,'Value'))


% --- Executes on button press in normfit.
function normfit_Callback(hObject, eventdata, handles)
% hObject    handle to normfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normfit
if ~isfield(handles,'data_normed')
    [fit,path] = uigetfile('*','Select fits');
    load([path fit]);
    
    acquisition=scd_data_separate(handles.data,handles.scheme);
    for i_acq=1:max(acquisition)
        handles.data_normed(:,:,:,acquisition==i_acq)=handles.data(:,:,:,acquisition==i_acq)./repmat(squeeze(fitted(:,:,:,1,6+i_acq)),[1 1 1 length(find(acquisition==i_acq))]);
    end
    handles.data_raw=handles.data;
end

if get(handles.normfit,'Value')
    handles.data=handles.data_normed;
else
    handles.data=handles.data_raw;
end

plotbutton_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes on button press in modelfitting.
function modelfitting_Callback(hObject, eventdata, handles)
% hObject    handle to modelfitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ax.scheme=handles.scheme(logical(handles.Selection),:);
Ax.Dr=str2double(get(handles.Dr,'String'));
Ax.Dcsf=str2double(get(handles.Dcsf,'String'));

if ~get(handles.csf,'Value'), Ax.Dcsf=0; end

handles.x =[];
% norm
Ax.data=double(squeeze(handles.data(handles.Y(1),handles.X(1),handles.Z,logical(handles.Selection))));
switch get(get(handles.panel_noise,'SelectedObject'),'Tag')
    case 'fixSNR'
        Ax.sigma_noise=1/str2double(get(handles.SNR,'String'))*max(Ax.data);
    case 'fixSigma'
        Ax.sigma_noise=str2double(get(handles.sigmanoise,'String'));
        set(handles.SNR,'String',num2str(max(Ax.data)/Ax.sigma_noise));
    case 'noisepervoxel'
        Ax.noisepervoxel=1;
end

models = get(handles.modelname,'String'); models = models{get(handles.modelname,'value')};
if strcmp(models,'CHARMED')
    Ax.optimizationfun = @scd_optimization_rician_likelihood;
    if get(handles.norm_fit,'Value'), Ax.norm='fit'; else Ax.norm='none'; end
    if get(handles.onediam,'Value'), Ax.onediam=1; else Ax.onediam=-1; end
    if get(handles.norm_fitT2,'Value'), Ax.fitT2 = 1; end
else
    Ax.optimizationfun = str2func(strrep(models,'.m',''));
    % create options
    opts = handles.opts;
    N=length(opts)/2;
    for i=1:N
        if islogical(opts{2*i})
            optionvalue = get(handles.modeloption_CUSTOM_handle(i),'Value');
        elseif isnumeric(opts{2*i})
            optionvalue = str2num(get(handles.modeloption_CUSTOM_handle(i),'String'));
        elseif iscell(opts{2*i})
            optionvalue = opts{2*i}{get(handles.modeloption_CUSTOM_handle(i),'Value')};
        end
        Ax.(matlab.lang.makeValidName(handles.opts{2*i-1}))=optionvalue;
    end
    
end

handles.Ax = Ax;

disp('<strong> Fitting voxels. please wait...</strong>')

for i_point=1:length(handles.Y)
    Ax.data=squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,logical(handles.Selection)));
    [handles.x(i_point,:), handles.modelfit{i_point}(logical(handles.Selection)),Ax_out] = Ax.optimizationfun(Ax);
end
for ip=1:length(Ax_out.parametersnames)
    disp(Ax_out.parametersnames{ip})
    disp(handles.x(:,ip)')
end

set(handles.GenerateMap,'enable','on')
set(handles.sigmanoise,'String',num2str(Ax_out.sigma_noise));


plotbutton_Callback(hObject, eventdata, handles)
handles.modelfit={};
guidata(hObject, handles);

% --- Executes on button press in norm_b0.
function norm_b0_Callback(hObject, eventdata, handles)
% hObject    handle to norm_b0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_b0


% --- Executes on button press in csf.
function csf_Callback(hObject, eventdata, handles)
% hObject    handle to csf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of csf


% --- Executes on button press in fixDh.
function fixDh_Callback(hObject, eventdata, handles)
% hObject    handle to fixDh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixDh


% --- Executes on button press in fitT2.
function fitT2_Callback(hObject, eventdata, handles)
% hObject    handle to fitT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fitT2


% --- Executes on button press in onediam.
function onediam_Callback(hObject, eventdata, handles)
% hObject    handle to onediam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onediam



function Dcsf_Callback(hObject, eventdata, handles)
% hObject    handle to Dcsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dcsf as text
%        str2double(get(hObject,'String')) returns contents of Dcsf as a double


% --- Executes during object creation, after setting all properties.
function Dcsf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dcsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Xaxis_qvalue.
function Xaxis_qvalue_Callback(hObject, eventdata, handles)
% hObject    handle to Xaxis_qvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Xaxis_qvalue
set(handles.Xaxis_bvalue,'Value',0)

% --- Executes on button press in Xaxis_bvalue.
function Xaxis_bvalue_Callback(hObject, eventdata, handles)
% hObject    handle to Xaxis_bvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Xaxis_bvalue
set(handles.Xaxis_qvalue,'Value',0)



function Dr_Callback(hObject, eventdata, handles)
% hObject    handle to Dr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dr as text
%        str2double(get(hObject,'String')) returns contents of Dr as a double
set(handles.Dr,'String',num2str(str2double(get(handles.Dr,'String'))));

% --- Executes during object creation, after setting all properties.
function Dr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SNR_Callback(hObject, eventdata, handles)
% hObject    handle to SNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR as text
%        str2double(get(hObject,'String')) returns contents of SNR as a double
set(handles.Dr,'String',num2str(str2double(get(handles.Dr,'String'))));


% --- Executes during object creation, after setting all properties.
function SNR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function Gmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Gmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gmax_edit as text
%        str2double(get(hObject,'String')) returns contents of Gmax_edit as a double


% --- Executes during object creation, after setting all properties.
function Gmax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in norm_fit.
function norm_fit_Callback(hObject, eventdata, handles)
% hObject    handle to norm_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_fit
if ~get(hObject,'Value'), set(handles.norm_fitT2,'Value',0); end

% --- Executes on button press in norm_fitT2.
function norm_fitT2_Callback(hObject, eventdata, handles)
% hObject    handle to norm_fitT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_fitT2
if get(hObject,'Value'), set(handles.norm_fit,'Value',1); end


% --- Executes on button press in gammadiam.
function gammadiam_Callback(hObject, eventdata, handles)
% hObject    handle to gammadiam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gammadiam



function sigmanoise_Callback(hObject, eventdata, handles)
% hObject    handle to sigmanoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmanoise as text
%        str2double(get(hObject,'String')) returns contents of sigmanoise as a double


% --- Executes during object creation, after setting all properties.
function sigmanoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmanoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotting_normalize.
function plotting_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to plotting_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotting_normalize


% --- Executes when selected object is changed in panel_noise.
function panel_noise_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_noise 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in discard.
function discard_Callback(hObject, eventdata, handles)
% hObject    handle to discard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
axes(handles.plot1)
selecttoolindex=selectdata_v2;
selectotherindex = find(handles.Selection);

handles.selecttoolindex = false(size(handles.Selection));
handles.selecttoolindex(selectotherindex(selecttoolindex)) = true;
guidata(hObject, handles);


% --- Executes on button press in GenerateMap.
function GenerateMap_Callback(hObject, eventdata, handles)
% hObject    handle to GenerateMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mkdir('qDLab_Maps')
Ax=handles.Ax; Ax=rmfield(Ax,'data'); Ax=rmfield(Ax,'scheme'); Select = ~~handles.Selection;
save(['qDLab_Maps' filesep 'fitting_param.mat'],'Ax')
save(['qDLab_Maps' filesep 'Selected_data'],'Select')
choice = questdlg('Create an output folder qDLab_Maps in the current folder?','Output directory','Yes','Choose folder','Yes');
if strcmp(choice,'Choose folder')
    outputfolder = uigetdir('.','Select output folder');
elseif isempty(choice)
    return;
else
    outputfolder = 'qDLab_Maps';
end
scd_diameter_map_parallel_computing(handles.data_fname_all,handles.scheme_fname_all,'output',outputfolder,'fitting_param','qDLab_Maps/fitting_param.mat','parallel',num2str(get(handles.UseParallel,'Value')),'mask',handles.mask_fname_all,'Select',Select)


% --- Executes on button press in UseParallel.
function UseParallel_Callback(hObject, eventdata, handles)
% hObject    handle to UseParallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseParallel


% --- Executes on button press in maskfname_button.
function maskfname_button_Callback(hObject, eventdata, handles)
% hObject    handle to maskfname_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[mask_fname,path] = uigetfile({'*.nii';'*.nii.gz'},'Select Mask'); handles.mask_fname_all=[path,mask_fname];
dims=size(handles.data(:,:,:,1));
if mask_fname
    nii = load_untouch_nii(handles.mask_fname_all);
    if min(size(nii.img)==dims)
        set(handles.maskfname,'String',mask_fname);
        guidata(hObject, handles);
    else
        disp('<strong>mask dimensions are different than your dataset</strong>')
    end
end


% --- Executes on selection change in modelname.
function handles = modelname_Callback(hObject, eventdata, handles)
% hObject    handle to modelname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modelname contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modelname
contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'CHARMED'
        set(handles.modeloption_CUSTOM,'visible','off');
        set(handles.modeloption_CHARMED,'visible','on');
    otherwise
        set(handles.modeloption_CHARMED,'visible','off');
        set(handles.modeloption_CUSTOM,'visible','on');
        
        % delete previous
        childs=get(handles.modeloption_CUSTOM,'Children');
        for ic=childs, delete(ic); end
        
        % Create buttons
        model=str2func(strrep(contents{get(hObject,'Value')},'.m',''));
        try opts=model(); 
        catch
            opts=[]; 
        end
        handles.opts = opts;
        if ~isempty(opts)
            N = length(opts)/2;
            [I,J]=ind2sub([4 3],1:2*N); Iw = 1/max(I); I=(I-1)/max(I); Jh = 1/max(J); J=(J-1)/max(J); J=1-J-Jh;
            for i = 1:N
                if islogical(opts{2*i})
                    hInput(i) = uicontrol('Style','checkbox','String',opts{2*i-1},...
                        'Parent',handles.modeloption_CUSTOM,'Units','normalized','Position',[[I(2*i)] J(2*i-1) Iw Jh],...
                        'Value',opts{2*i},'HorizontalAlignment','center');
                elseif isnumeric(opts{2*i})
                    uicontrol('Style','Text','String',opts{2*i-1},...
                        'Parent',handles.modeloption_CUSTOM,'Units','normalized','Position',[I(2*i-1) J(2*i-1) Iw Jh]);
                    hInput(i) = uicontrol('Style','edit',...
                        'Parent',handles.modeloption_CUSTOM,'Units','normalized','Position',[I(2*i) J(2*i) Iw Jh],'String',opts{2*i});
                elseif iscell(opts{2*i})
                    uicontrol('Style','Text','String',opts{2*i-1},...
                        'Parent',handles.modeloption_CUSTOM,'Units','normalized','Position',[I(2*i-1) J(2*i-1) Iw Jh]);
                    hInput(i) = uicontrol('Style','popupmenu',...
                        'Parent',handles.modeloption_CUSTOM,'Units','normalized','Position',[I(2*i) J(2*i) Iw Jh],'String',opts{2*i});
                    
                end
            end
            handles.modeloption_CUSTOM_handle=hInput;
        end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function modelname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in XaxisCustom.
function XaxisCustom_Callback(hObject, eventdata, handles)
% hObject    handle to XaxisCustom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XaxisCustom


% --- Executes during object creation, after setting all properties.
function text33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when option_CHARMED is resized.
function option_CHARMED_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to option_CHARMED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when modeloption_CHARMED is resized.
function modeloption_CHARMED_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to modeloption_CHARMED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on selection change in noddimodellist.
function noddimodellist_Callback(hObject, eventdata, handles)
% hObject    handle to noddimodellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns noddimodellist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from noddimodellist



% --- Executes during object creation, after setting all properties.
function noddimodellist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noddimodellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% modellist=matlab.codetools.requiredFilesAndProducts('SynthMeas.m','toponly');
% modellist(1)=[];
% for im=1:length(modellist)
%     [~,modellist{im}] = fileparts(modellist{im});
%     modellist{im} = strrep(modellist{im},'SynthMeas','');
% end
% set(hObject,'String',modellist)
