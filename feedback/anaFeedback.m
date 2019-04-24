function [out, file]=anaFeedback(file)
% Current only do this for dBz scans
global fbdata
if ~exist('file','var'), file = getFiles('*dBz*'); file = fliplr(file); end
out = anadBz(file);
scanDate=[];
for i = 1:length(out.a)
    out.nloop(i) = out.a(i).scan.data.conf.nloop;
    out.pulseTime = out.a(i).scan.data.pulsegroups(1).params(1) * out.a(i).scan.data.pulsegroups(1).npulse(1);
    if isfield(out.a(i).scan.data,'gradHist')
        out.setPt(i) = out.a(i).scan.data.setpt(1);        
        gradHist{i} = [out.a(i).scan.data.gradHist{:}];        
        meanGradData(:,i) = nanmean(gradHist{i},2);
        stdGradData(:,i) = nanstd(gradHist{i},[],2);
        out.nPump(i) = sum(out.a(i).scan.data.FBset);
        gradErr{i} = gradHist{i}(1,:)-out.setPt(i);
        scanDate(i) =  out.scanDate(i); 
    else
        meanGradData(1:3,i)=nan; 
        stdGradData(1:3,i)=nan; 
        out.nPump(i) = nan;         
    end
end
figure(14); clf;
pumptime=fbdata.pumptime(1)*1e3; % units of ms,% In the future, let's save this with data. 
ha = tightSubplot([2,2]);

plot(ha(3),out.nloop.*out.pulseTime/1e3,out.t2s,'.');
xlabel(ha(3),'Measurement Time (ms)'); ylabel(ha(3),'T_2*'); 

plot(ha(1),meanGradData(2,:)/pumptime,'.','DisplayName','Singlet'); hold(ha(1),'on');
plot(ha(1),meanGradData(3,:)/pumptime,'.','DisplayName','Triplet');
ylabel(ha(1),'Pump Rate (MHz/ms)'); legend(ha(1),'show'); 
plot(ha(2),stdGradData(2,:)/pumptime','.','DisplayName','Singlet'); hold(ha(2),'on');
plot(ha(2),stdGradData(3,:)/pumptime','.','DisplayName','Triplet');
ylabel(ha(2),'std Pump Rate (MHz/ms)'); legend(ha(2),'show'); 

plot(ha(4),scanDate,'.-'); 
end