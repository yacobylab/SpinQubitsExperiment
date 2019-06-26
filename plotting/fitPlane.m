function coeff = fitPlane(data)
[gx,gy] = gradient(data);
for i = 1:size(gx,1)
    gx(i,:)=smooth(gx(i,:),3);
end
for i = 1:size(gy,2)
    gy(:,i)=smooth(gy(:,i),3);
end
coeff(1)=nanmedian(cull(gx));
coeff(2)=nanmedian(cull(gy));
coeff(3)=nanmean(nanmean(data));
end

function filtData=cull(data)
med = nanmedian(data(:));
medDist = abs(data(:)-med); 
medSD = nanmedian(medDist);
filtData = data(medDist < 2 * medSD);
meanData = nanmean(filtData);
filtData = data(abs(data(:)-meanData) < 2*medSD);
end