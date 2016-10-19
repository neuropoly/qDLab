function Ib0 = scd_preproc_getIb0(data,scheme)
TEs=unique(scheme(:,7));
Ib0=zeros(size(data));
for i=1:length(TEs)
    indexTEi=scheme(:,7)==TEs(i);
    bmin = min(scd_scheme2bvecsbvals(scheme(indexTEi,:)));
    if bmin>300, warning(['no b=0 data for TE = ' num2str(TEs(i))]), end
    if sum(size(data)>1)==1
        dataTE=data(indexTEi);
        Ib0(indexTEi)=mean(dataTE(scd_scheme2bvecsbvals(scheme(indexTEi,:))==bmin));
    else
        dataTE=data(:,:,:,indexTEi);
        Ib0(:,:,:,indexTEi)=repmat(mean(dataTE(:,:,:,scd_scheme2bvecsbvals(scheme(indexTEi,:))==bmin),4),[1 1 1 sum(indexTEi)]);
    end
end