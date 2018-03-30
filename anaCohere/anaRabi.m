function [ yvals,rabi,d ] = anaRabi(fname)
% analyzes rabi frequency vs drive with the same pulsegroup
%[yvals,rabi,d ] = anaRabi(fname)
%   fname: optional filename

if ~exist('fname','var'), fname=get_files('sm*.mat'); end
d=ana_avg(fname,struct('ops',''));

data=(d.data{1});
dataAvg=squeeze(nanmean(data));
try
    yvals=linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);
catch
    yvals=d.tv;
end
xvals=d.xv{1};

%try to see if its a multi file
ind=find(diff(xvals)<0);
if ~isempty(ind)
    xvals=xvals(ind+1:end);
    dataAvg=dataAvg(:,ind+1:end);
end
figure(1); clf;
subplot(4,1,1); imagesc(xvals,yvals,dataAvg); set(gca,'YDir','norm');
xlabel('Burst time (ns)'); ylabel('Rabi amplitude'); colorbar;

% plot and fit a slice
figure(200); clf;
sliceStart=1;        
for i=1:size(dataAvg,1)
    yvals=dataAvg(i,:);    
    fp = fioscill(xvals, yvals, 2); %estimate the freq and phase    
    xvals=xvals(sliceStart:end);
    yvals=yvals(sliceStart:end);   
        
    fitfn=@(p,x) p(1)+p(2)*cos(2*pi*x.*(p(3))+p(4)).*exp(-(x./p(5)).^p(6))+p(7).*x;
    beta=[mean(yvals) (max(yvals)-min(yvals))/2 fp(4)/(2*pi) 0 xvals(end)/2 2 0];    
    mask=[1 1 1 1 1 0 0];
    [beta,~,~,~,~,err]=fitwrap('plinit plfit',xvals,yvals,beta,fitfn,mask);    
    fit=fitfn(beta,xvals);    
    figure(200);  hold on; plot(xvals,yvals+i); plot(xvals,fit+i);    
    rabi(i)=abs(beta(3))*1e3; T2(i)=abs(beta(5)); T2err(i)=err(1,5,1); amp(i)=beta(2);
end

figure(1); subplot(4,1,2);
errorbar(yvals,T2,T2err,'-o'); xlabel('Rabi amplitude'); ylabel('T2*');
subplot(4,1,3);
plot(yvals,rabi,'o-');
xlabel('Rabi amplitude'); ylabel('Rabi Freq');
subplot(4,1,4);
plot(yvals,amp,'o-');
xlabel('Rabi amplitude'); ylabel('Visiblity');
end