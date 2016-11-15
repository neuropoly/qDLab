function scd_display_fits(data,model,scheme)
% scd_display_fits(data,model,scheme)
% Plot model vs data curves. 
%
% Inputs:
% data          vector(nb_dwi)          Experimental Diffusion Mri Signal
% model         vector(nb_dwi)          Model Diffusion Mri Signal
% scheme        vector(nb_dwi,9)        Diffusion Mri Parameters (loaded with scd_schemefile_read)
%  
%
% Example:
% img=load_nii_data('diffusion.nii');
% data=img(x,y,z,:)
% scheme=scd_schemefile_read('acq.scheme');
% data_model = x(8)*exp(-scheme(:,7)./x(9)).*scd_model_GPD_composite(x,Ax);
% scd_display_fits(data,data_model,scheme)
%
% See Also: scd_display_qspacedata

scd_display_qspacedata(data,scheme)
hold on
scd_display_qspacedata(model,scheme,0,'none','-')
ylim([0 max(max(data), max(model))])

