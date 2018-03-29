function [flatData,slp]=flattenData(data,xvalsM) 
% flatten the data on each line using 4th order poly fit. 
% function [flatData,slp]=flattendata(data,xvalsM) 
szData = size(data); 
if ~exist('xvalsM','var') 
xvalsM = 1:szData(2); 
end
xChar = [xvalsM(10),mean(xvalsM),xvalsM(end-10)]; 
wid = 'MATLAB:polyfit:PolyNotUnique'; warning('off',wid);
for i = 1:szData(1)  
    xvals = xvalsM; 
    yvals = data(i,:); 
    inds = find(isnan(yvals)); 
    xvals(inds)=[]; yvals(inds)=[];
    poly = 1; 
    
        if poly 
            p = polyfit(xvals, yvals, 4);
            flatData(i,:) = data(i,:)-p(1)*xvalsM.^4-p(2)*xvalsM.^3-p(3)*xvalsM.^2-p(4)*xvalsM-p(5);                 
            slp(i,:) = 3*p(1)*xChar.^2 + 2 * p(2)*xChar + p(3); 
        else
            
            p = polyfit(xvals, yvals, 1);
            flatData(i,:) = data(i,:)-p(1)*xvalsM-p(2);  
            
            %sinFit = @(p,x) p(1).*sin(p(2).*x+p(3))+p(4).*x+p(5); 
%             amp = std(flatData(i,:)); phase = 0.1; omega = range(xvalsM)*2*pi*4; 
            %beta0 = [amp,omega,phase,p(1),p(2)];     
%             sinFit = @(p,x) p(1).*sin(p(2).*x+p(3)); 
%             beta0=[amp,omega,phase]; 
%             params=fitwrap('plinit plfit',xvalsM,flatData(i,:),beta0,sinFit); 
            %params=fitwrap('plinit plfit',xvalsM,data(i,:),beta0,sinFit); 
            %flatData(i,:)=data(i,:)-sinFit(params,xvalsM); 
            %flatData(i,:)=flatData(i,:)-sinFit(params,xvalsM); 
            slp(i,:) = p(1); 
        end        
end


end