function [out,file] = anadBz(file)
% Process, fit and plot dBz data.
% function anadBz
% Fit to sin decaying as t^2
% print T2, return list of parameters: offset, out.amplitude, out.freq (rad/s), out.phase, t2*

if ~exist('file','var')
    file = getFiles('*dBz*'); file = fliplr(file); 
end
a = procPlsData(file,'noplot nodbz');
off = 0.25;
% offset, out.amplitude, out.freq (rad/s), out.phase, t2*)

%a = fliplr(a);
out.a = a;
out.scanDate = datetime([a.scantime],'ConvertFrom','datenum'); 
figure(17); clf; hold on;
for i = 1:length(a)    
    data = squeeze(nanmean(a(i).data{1}))+off*(i-1);    
    [pars,~,~,res,mse,err] = fitosc(a(i).xv{1},data','pldata samecolor phase');
    out.t2s(i) = 1./pars(6); out.t2sErr(i) = out.t2s(i)^2*err(6);
    out.amp(i) = pars(2); out.phase(i) = pars(3);
    out.ampErr(i) = err(2); out.phaseErr(i) = err(3);
    out.freq(i) = pars(4)/2/pi*1e3; out.freqErr(i) = err(4)/2/pi*1e3;
    fid(i) = a(i).fidelity; 
    if length(a)==1
        fprintf('T2*: %3.1f +/- %3.1f ns. Amplitude %1.2f \n',out.t2s,out.t2sErr,out.amp);
        fprintf('Frequency: %3.3f +/- %3.3f MHz \n',out.freq,out.freqErr);
        fprintf('T1 = %3.3g us. Peak separation %3.3g mV. Fidelity %3.2f%% \n',a(i).t1*1e6,1e3*diff(a(i).meanvals),a(i).fidelity*100);
    end
end

if length(a) > 1
    figure(99); clf;
    ha=tightSubplot([2,3]);
    errorbar(ha(1),out.t2s,out.t2sErr);
    ylabel(ha(1),'T_2* (ns)');
    errorbar(ha(2),out.phase/pi,out.phaseErr);
    ylabel(ha(2),'Phase');
    errorbar(ha(3),out.amp,out.ampErr);
    ylabel(ha(3),'Amplitude');
    errorbar(ha(4),out.freq,out.freqErr);
    ylabel(ha(4),'Frequency (MHz)');
    plot(ha(5),fid*100,'.-');
    ylabel(ha(5),'Fidelity');
    plot(ha(6),out.scanDate,'.-')
end
end