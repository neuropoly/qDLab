function [] = scd_display_qspacedata(data,scheme,Marker,Linestyle,color,noise)
% scd_display_qspacedata(data,scheme)
%data: vector OR NxM matrix. N: number of measurement (=size(scheme,1)). M: number of
%      observation (for error bars).

%  selection=ismember(scheme(:,9),selection);
% scheme_selected=scheme(selection,:);

if ~exist('Marker','var'), Marker='x'; end
if ~exist('Linestyle','var'), Linestyle='none'; end

q=double(scheme(:,8));
[q,I] = sort(q);
data = data(I); scheme = scheme(I,:);


seq=unique(scheme(:,9)); ND=length(seq);
if ~exist('color','var')
    color=jet(ND);
end

for iD=1:ND
    seqiD=find(scheme(:,9)==seq(iD));
    
    datavoxel=mean(squeeze(data),2);
    
    
    g(iD)=plot(q(seqiD),datavoxel(seqiD),'LineStyle',Linestyle, 'Marker',Marker,'Color',color(min(iD,end),:),'LineWidth',2);
    hold on
    
    set(g(iD),'DisplayName',['Delta=' num2str(mean(scheme(seqiD,5))) ' delta=' num2str(mean(scheme(seqiD,6))) ' TE=' num2str(mean(scheme(seqiD,7)))]);
    
    
    if exist('noise','var')==1
        
        h(iD)=errorbar(q(seqiD),datavoxel(seqiD),noise(seqiD), 'xr', 'Color', color);
    elseif size(data,2)>1
        data_stdvoxel=std(data,0,2);
        h(iD)=errorbar(q(seqiD),datavoxel(seqiD),data_stdvoxel(seqiD), 'xr', 'Color', color);
        %         figure(3)
        %         d(iD)=plot(q(seqiD),data_stdvoxel(seqiD),'LineStyle','none', 'Marker','x','Color',color,'LineWidth',2);
        %         hold on
    end
    if exist('h','var')==1
        % don't display data legend
        hAnnotation = get(h(iD),'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off');
    end
end


xlabel('q','FontSize',15); ylabel('Signal (%b0)','FontSize',15);


legend('show')
set(gca,'FontSize',15)
%ylim([0 1.2])
grid on, box off

hold off

end
