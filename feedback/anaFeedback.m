function [out, file]=anaFeedback(file)
% Current only do this for dBz scans
if ~exist('file','var'), file = getFiles('*dBz*'); file = fliplr(file); end
out = anadBz(file);
scanDate=[];
for i = 1:length(out.a)
    out.nloop(i) = out.a(i).scan.data.conf.nloop;
    out.pulseTime = out.a(i).scan.data.pulsegroups(1).params(1) * out.a(i).scan.data.pulsegroups(1).npulse(1);
    if isfield(out.a(i).scan.data,'gradHist') && ~isempty(out.a(i).scan.data.gradHist)
        out.setPt(i) = out.a(i).scan.data.setpt(1);        
        gradHist{i} = [out.a(i).scan.data.gradHist{:}];
        pumpHist{i} = [out.a(i).scan.data.pumpHist{:}];
        meanGradData(:,i) = nanmean(gradHist{i},2);
        stdGradData(:,i) = nanstd(gradHist{i},[],2);
        gradErr{i} = gradHist{i}(1,:)-out.setPt(i);        
        out.nPump(i) = sum(out.a(i).scan.data.FBset);        
        scanDate{i} =  out.scanDate(i); 
        ind = str2double(out.a(i).scan.loops(1).getchan{1}(end)); 
        try
            pumptime(i) = out.a(i).scan.data.FBData.params(ind).pumptime*1e3; 
        catch
            pumptime(i) = out.a(i).scan.data.FBData.pumptime(ind)*1e3; 
        end
    else
        out.t2s(i) = nan; 
        meanGradData(1:3,i)=nan; 
        stdGradData(1:3,i)=nan; 
        out.nPump(i) = nan;         
    end
end
figure(14); clf;

ha = tightSubplot([2,2]);

plot(ha(3),out.nloop.*out.pulseTime/1e3,out.t2s,'.');
xlabel(ha(3),'Measurement Time (ms)'); ylabel(ha(3),'T_2*'); 

plot(ha(1),abs(meanGradData(2,:)./pumptime),'.','DisplayName','Singlet'); hold(ha(1),'on');
plot(ha(1),abs(meanGradData(3,:)./pumptime),'.','DisplayName','Triplet');
ylabel(ha(1),'Pump Rate (MHz/ms)'); legend(ha(1),'show'); 
plot(ha(2),stdGradData(2,:)./pumptime','.','DisplayName','Singlet'); hold(ha(2),'on');
plot(ha(2),stdGradData(3,:)./pumptime','.','DisplayName','Triplet');
ylabel(ha(2),'std Pump Rate (MHz/ms)'); legend(ha(2),'show'); 

%plot(ha(4),scanDate,'.-'); 
formatFig(14,'exch full',2,2)
end