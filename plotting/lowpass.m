function data=lowpass(sigma,data)
% Spyview style lowpass gaussian filter. 
% function data=lowpass(sigma,data)
% sigma is standard deviation of indices to average over. 
if sigma(1) == 0 && sigma(2) == 0
    return;
end
xs=ceil(3*sigma(1))*2+1;
ys=ceil(3*sigma(2))*2+1;
if sigma(1) == 0
    y=((-floor(ys/2)):1:floor(ys/2))';
    kernel=exp(-y.^2/(2*sigma(2)));
elseif sigma(2) == 0
    x=(-floor(xs/2)):1:floor(xs/2);
    kernel=exp(-x.^2/(2*sigma(1)));
else
    y=((-floor(ys/2)):1:floor(ys/2));
    x=(-floor(xs/2)):1:floor(xs/2);
    [cx,cy]=meshgrid(x,y);
    kernel=exp(-cx.^2/(2*sigma(1)) + -cy.^2/(2*sigma(2)));
end
kernel = kernel / sum(kernel(:));
pmp=padarray(data,[floor(ys/2),floor(xs/2)],'replicate');
pmp=conv2(pmp,kernel,'same');
data=pmp(1+floor(ys/2):(end-floor(ys/2)),1+floor(xs/2):(end-floor(xs/2)));
end
