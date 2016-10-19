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

% Last Modified by GUIDE v2.5 19-Oct-2016 17:31:02

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

handles.selecttoolindex = false([size(handles.scheme,1) 1]);

% Simplify options depending on the dataset:
    % fit T2?
if length(unique(handles.scheme(scd_scheme2bvecsbvals(handles.scheme)<1000,7))) == 1 % Don't propose to fit T2 if only one echo time
    set(handles.norm_fitT2,'Visible','off');
else
    set(handles.norm_fitT2,'Value',true);
end

    % fit diameter distribution?
if ~sum(scd_scheme2bvecsbvals(handles.scheme)>30000), set(handles.gammadiam,'enable','off'); end % don't propose gamma distribution if maximal bvalue<30,000

    % estimate noise voxel-wise?
[~,c]=consolidator(handles.scheme(:,1:8),[],'count');
cmax = max(c); % find images repeated more than 5 times (for relevant STD)
if cmax<5
    set(handles.noisepervoxel,'enable','off')
elseif cmax>15 % noise using STD is robust --> select this option by default
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
varargout{1} = get(handles.AcquisitionList,'UserData');


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
    [scheme_fname,path] = uigetfile('*','Select schemefile');scheme_fname_all=[path,scheme_fname];
else
    scheme_fname = scheme_fname_all;
end
set(handles.textscheme,'String',['SchemeFile : ' scheme_fname])
handles.scheme_fname_all=scheme_fname_all;

if scheme_fname
    handles.scheme=scd_schemefile_read(scheme_fname_all);
    if size(handles.data,4)~=size(handles.scheme,1), error(['<strong>Error: your dataset has ' num2str(size(handles.data,4)) ' while your schemefile has ' num2str(size(handles.scheme,1)) ' rows.</strong>']); end

    % sort data
    [handles.q,Isort]=sort(handles.scheme(:,8));
    handles.data=handles.data(:,:,:,Isort); handles.scheme=handles.scheme(Isort,:);
    % list acquisitions
    handles.acqList=unique(handles.scheme(:,7:-1:5),'rows');
    for i=1:size(handles.acqList,1), acqListString{i}=sprintf('% 15g % 15g % 15g',handles.acqList(i,:)); end
    set(handles.AcquisitionList,'String',acqListString);
    set(handles.AcquisitionList,'Value',1:size(handles.acqList,1)); % select all by default
    % list directions
    handles.bvecs=unique(handles.scheme(:,1:3),'rows');
    for i=1:size(handles.bvecs,1), dirList{i}=sprintf('% 15g % 15g % 15g',handles.bvecs(i,:)); end
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
    [data_fname,path] = uigetfile({'*.nii';'*.nii.gz'},'Select data');handles.data_fname_all=[path,data_fname];
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
    handles.Selection=handles.Selection + i.*(dir_logical & handles.scheme(:,7)==handles.acqList(acqselected(i),1) & handles.scheme(:,6)==handles.acqList(acqselected(i),2) & handles.scheme(:,5)==handles.acqList(acqselected(i),3) & ~handles.selecttoolindex);
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
    Ax=scd_scheme2struct(handles.scheme);
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
        [qvalues,ia,iq]=unique(schemeiaq(:,8)); schemeiaq = schemeiaq(ia,:);
            % abscissa
        if get(handles.bvalue_radio,'Value'), bvals=scd_scheme2bvecsbvals(handles.scheme(handles.Selection==iaq,:)); absc=bvals(ia); else absc=qvalues; end
        [~,data]=consolidator(iq,squeeze(handles.data(handles.Y(i_point),handles.X(i_point),handles.Z,handles.Selection==iaq)));
            % normalize data
            if get(handles.plotting_normalize,'value')
                if isempty(mdel)
                    normvalue=scd_preproc_getIb0(data,schemeiaq);
                else
                    normvalue = max(mdel(handles.Selection==iaq));
                end
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
            set(handles.h(plotnumber),'DisplayName',['Delta=' num2str(handles.acqList(acqselected(iaq),3)) ' delta=' num2str(handles.acqList(acqselected(iaq),2)) ' TE=' num2str(handles.acqList(acqselected(iaq),1))]);
        else
            hAnnotation = get(handles.h(plotnumber),'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
        end
        
        
        % plot model
        if ~isempty(mdel)
                 % abscissa
                if get(handles.bvalue_radio,'Value'), bvals=scd_scheme2bvecsbvals(handles.scheme(handles.Selection==iaq,:)); absc=bvals; else absc=handles.q(handles.Selection==iaq); end

                handles.g(iaq)=plot(absc,mdel(handles.Selection==iaq)./normvalue,'Color',handles.colorplot(plotnumber,:));
                set(handles.g(iaq),'Linewidth',3)       
                hAnnotation = get(handles.g(iaq),'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')

        end


    end
end
ylim('auto')
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
function Zslider_Callback(hObject, eventdata, handles)
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
if get(handles.norm_b0,'Value'), Ax.norm='fit'; else Ax.norm='none'; end
if get(handles.csf,'Value'), Ax.Dcsf=str2double(get(handles.Dcsf,'String')); end
if get(handles.onediam,'Value'), Ax.onediam=1; else Ax.onediam=0; end
Ax.Dr=str2double(get(handles.Dr,'String'));

switch get(get(handles.uipanel4,'SelectedObject'),'Tag')
    case 'norm_b0'
        Ax.norm='none';
    case 'norm_fit'
        Ax.norm='fit';
    case 'norm_fitT2'
        Ax.norm='none';
        Ax.fitT2=1;
end
        

handles.x =[];
% norm
Ax.data=squeeze(handles.data(handles.Y(1),handles.X(1),handles.Z,logical(handles.Selection)));
switch get(get(handles.panel_noise,'SelectedObject'),'Tag')
    case 'fixSNR'
        Ax.sigma_noise=1/str2double(get(handles.SNR,'String'))*max(Ax.data);
    case 'fixSigma'
        Ax.sigma_noise=str2double(get(handles.sigmanoise,'String'));
        set(handles.SNR,'String',num2str(max(Ax.data)/Ax.sigma_noise));
    case 'noisepervoxel'
        Ax.noisepervoxel=1;
end


if get(handles.modelname,'value') == 1, 
    Ax.optimizationfun = @scd_optimization_rician_likelihood;
elseif get(handles.modelname,'value') == 2,
    Ax.optimizationfun = @scd_optimization_custom_model;
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


% --- Executes on button press in qvalue_radio.
function qvalue_radio_Callback(hObject, eventdata, handles)
% hObject    handle to qvalue_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of qvalue_radio
set(handles.bvalue_radio,'Value',0)

% --- Executes on button press in bvalue_radio.
function bvalue_radio_Callback(hObject, eventdata, handles)
% hObject    handle to bvalue_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bvalue_radio
set(handles.qvalue_radio,'Value',0)



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


% --- Executes on button press in norm_fitT2.
function norm_fitT2_Callback(hObject, eventdata, handles)
% hObject    handle to norm_fitT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_fitT2


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
mkdir('AxCaliber_Maps')
Ax=handles.Ax; Ax=rmfield(Ax,'data'); Ax=rmfield(Ax,'scheme'); Select = ~~handles.Selection;
save AxCaliber_Maps/fitting_param.mat Ax
save AxCaliber_Maps/Selected_data Select
scd_diameter_map_parallel_computing(handles.data_fname_all,handles.scheme_fname_all,'output','AxCaliber_Maps','fitting_param','AxCaliber_Maps/fitting_param.mat','parallel',num2str(get(handles.UseParallel,'Value')),'mask',handles.mask_fname_all,'Select',Select)


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
set(handles.maskfname,'String',mask_fname);
guidata(hObject, handles);


% --- Executes on selection change in modelname.
function modelname_Callback(hObject, eventdata, handles)
% hObject    handle to modelname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modelname contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modelname


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
