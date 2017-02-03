function [] = scd_display_qspacedata_3D(data, scheme)
%chargement des données nifti


Angle=double(abs(scheme(:,3)));


% get point


h=unique(scheme(:,9));
ND=length(h);






if ~exist('color','var')
    color=jet(ND);
end





for i=1:ND
    
    
    
  if  size(data,2)==1 
    data=squeeze(data);
    figure(2)
    g(i)=plot(Angle(scheme(:,9)==h(i)),Voxel(scheme(:,9)==h(i)), 'LineStyle','none', 'Marker','x', 'Color',color(i,:), 'MarkerSize',10);
    hold on
    
    
  
  else
      
        figure(2)
        
        datavoxel=mean(data, 2);

        data_stdvoxel=std(data,0,2);
        
        
        g(i)=plot(Angle(scheme(:,9)==h(i)),datavoxel(seqiD),'LineStyle','none', 'Marker','x','Color',color,'LineWidth',2);
        hold on
  
    
    % don't display data legend
        hAnnotation = get(g(iD),'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        
        g(i)=errorbar(q(seqiD),datavoxel(seqiD),data_stdvoxel(seqiD), 'xr', 'Color', color);
        
        figure(3)
        d(iD)=plot(Angle(scheme(:,9)==h(i)),data_stdvoxel(seqiD),'LineStyle','none', 'Marker','x','Color',color,'LineWidth',2);
        hold on
  
end





end


title('E=f(Gz/G)', 'FontSize', 14);
xlabel('Gz/G ', 'FontSize', 14); 
ylabel('E ', 'FontSize', 14); 




end
