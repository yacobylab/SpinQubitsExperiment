function [out, file]=anaFeedback2(file)
% Current only do this for dBz scans
if ~exist('file','var'), file = getFiles; file = fliplr(file); end
out.a = procPlsData(file);
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
        ind = str2double(out.a(i).scan.loops(1).getchan{1}(end)); 
        try
            pumptime(i) = out.a(i).scan.data.FBData.params(ind).pumptime*1e3; 
        catch
            pumptime(i) = out.a(i).scan.data.FBData.pumptime(ind)*1e3; 
        end
        figure(41); clf;
        a = tightSubplot([2,2]); 
        hold(a(1),'on'); hold(a(2),'on');
        for j = 2:length(out.a(i).scan.data.gradHist)
            if ~isempty(out.a(i).scan.data.gradHist{j})
                plot(a(2),out.a(i).scan.data.err{j},'.-')
                plot(a(1),out.a(i).scan.data.pumpHist{j}(1,:),'.-')
                fGrad(j) = out.a(i).scan.data.err{j}(1,end); 
            end
        end
        plot(a(3),fGrad,'.'); 
    else
        out.t2s(i) = nan; 
        meanGradData(1:3,i)=nan; 
        stdGradData(1:3,i)=nan; 
        out.nPump(i) = nan;         
    end
end
figure(14); clf;

ha = tightSubplot([2,2]);

plot(ha(1),abs(meanGradData(2,:)/pumptime),'.','DisplayName','Singlet'); hold(ha(1),'on');
plot(ha(1),abs(meanGradData(3,:)/pumptime),'.','DisplayName','Triplet');
ylabel(ha(1),'Pump Rate (MHz/ms)'); legend(ha(1),'show'); 
plot(ha(2),stdGradData(2,:)/pumptime','.','DisplayName','Singlet'); hold(ha(2),'on');
plot(ha(2),stdGradData(3,:)/pumptime','.','DisplayName','Triplet');
ylabel(ha(2),'std Pump Rate (MHz/ms)'); legend(ha(2),'show'); 

%plot(ha(4),scanDate,'.-'); 
formatFig(14,'exch full',2,2)
end