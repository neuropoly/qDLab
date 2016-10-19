function scheme = scd_schemefile_read(scheme_file,varargin)
% load scheme files
% scd_schemefile_read('acq.scheme',(units in SI?))
% output: g_x  g_y  g_z  |G|(mT/um) Delta(ms) delta(ms) TE(ms) q(µm-1) acq#
% fid = fopen(scheme_file,'r');
% scheme=cell2mat(textscan(fid,'%f %f %f %f %f %f %f','delimiter','\n','CommentStyle','#'));

scheme=txt2mat(scheme_file,'InfoLevel',0);
scheme = round(scheme*10^6)/10^6; % floating-point problems

if isempty(varargin)
    % convert units
    if sum(scheme(:,1:3))>0
        scheme(:,4)=scheme(:,4).*sqrt(sum(scheme(:,1:3).^2,2))*1e-3; % G mT/um
    else
        scheme(:,4)=scheme(:,4).*1e-3; % G mT/um
    end
    scheme(:,1:3)=scheme(:,1:3)./repmat(sqrt(scheme(:,1).^2+scheme(:,2).^2+scheme(:,3).^2),1,3); scheme(isnan(scheme))=0;
    scheme(:,5) = scheme(:,5)*10^3; % DELTA ms
    scheme(:,6) = scheme(:,6)*10^3; % delta ms
    scheme(:,7) = scheme(:,7)*10^3; % TE ms
    gyro = 42.57; % kHz/mT
    if size(scheme,2)==7
        scheme(:,8) = gyro*scheme(:,4).*scheme(:,6); % um-1
    end
    % Detect different sequences
    list=unique(scheme(:,7:-1:5),'rows');
    nnn = size(list,1);
    for j = 1 : nnn
        for i = 1 : size(scheme,1)
            if  scheme(i,7:-1:5) == list(j,:)
                scheme(i,9) = j;
            end
        end
    end
       
end

disp([num2str(size(scheme,1)) ' images; ' num2str(length(find(scheme(:,4)<1e-8))) ' b=0 images; ' num2str(size(unique((scheme(:,9))),1)) ' different combinations Delta/delta/TE'])