function result = scd_readresults(data_file, output, onediam, scheme)
% read results from scd_diameter_map_parallel_computing
% result = scd_readresults(data_file, output, onediam, scheme)
% Example: scd_readresults('qspace_crop_eddy_moco.nii', 'ResultsAxCaliber/', 1, scd_schemefile_read('336810.scheme'))

if ~isdeployed, dbstop if error; end

result=load([output filesep 'fitted_results.mat']);

fields=result.parametersnames;

for ifield=1:length(fields)
    if ~isempty(fields{ifield})
        if isstr(data_file) && ~isempty(data_file)
            save_nii_v2((mean(result.fitted(:,:,:,:,ifield),4)),[output filesep fields{ifield}],data_file,64);
        end
    end
end

% save fr
if isfield(fields,'fh')
    dims = size(result.fitted);
    fh   = (mean(result.fitted(:,:,:,:,strcmp(fields,'fh')),4));
    fcsf = (mean(result.fitted(:,:,:,:,strcmp(fields,'fcsf')),4));
    if isempty(fcsf), fcsf=zeros(dims(1:3)); end
    fh(~fh)=1;
    if isstr(data_file) && ~isempty(data_file)
        save_nii_v2(ones(dims(1:3)) - fh - fcsf,[output filesep 'fr'],data_file);
    end
end
end