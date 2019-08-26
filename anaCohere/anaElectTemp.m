%files=getFiles('sm_leadPos*');
d = loadFiles; 
files=fliplr(files);

MCtemp=[];
etemp=[]; arm=[]; data = [];

beta0=[.029 .15];
mask=[1 1];
%amp=10e-3;
width=1; % mV on either side of center to 

widthFit=[]; widthErr=[]; slopeFit=[]; centerFit=[];
offstFit=[]; ampFit=[]; unitlessTemp=[];

figure(5); clf; ha = tightSubplot(length(d),'smart title'); 
for i = 1:length(d)    
    data=d(i).data{1}(12,:);
    %data=squeeze(nanmean(d.data{1}));
    try
        MCtemp(i)=1e3*d(i).configvals(strcmp(d(i).configch,'MC')); %#ok<*SAGROW>
    catch
        %get the temperature from the title.
        ind=regexp(files{i},'mK');
        temp=files{i}(ind-3:ind-1);
        temp(temp=='_')=[];
        MCtemp(i)=str2double(temp)*1e-3;
    end
    % xv has units of mV.
    xv = d(i).scan.data.pulsegroups(1).varpar; 
    xv = xv(:,2); %sqrt(xv(:,1).^2 + xv(:,2).^2); 
    % make an educated case where the center is
    [~, ind] = max(diff(data));
    center=xv(ind);
    inds = xv > (center-width) & xv < (center+width);
    xv = xv(inds);
    data=data(inds);
    xv = flipud(xv);
    data = flipud(data);

    fitfn=@(p,x) p(1) +p(2).*x + p(3)./(exp(((x-p(4))./p(5)))+1)+p(6).*x.^2; %fermi
    %fitfn=@(p,x) p(1) +p(2).*x -p(3).*tanh((x-p(4))./p(5))+p(6).*x.^2; %tanh   
    
    mask = [1 1 1 1 1 1]; 
    slp = (data(2)-data(1))./(xv(2)-xv(1)); 
    %       tc, alpha, offset, slope, amp, center,     
    %beta0 = [0.01, 0.1, data(2)-4*slp, slp, slp*4, center, MCtemp(i), 0, 0]; 
    %      offset,      slp, amp,  center,    quad. back.
    beta0=[data(1)-slp*2 slp slp*4 center .1 0];
    
    [params,~,~,~,~,err] = fitwrap('plinit plfit',xv',data,beta0,fitfn);
    fit=fitfn(params,xv);
    plot(ha(i),xv,data,'k.'); hold(ha(i),'on'); 
    plot(ha(i),xv,fit,'b');
    
    widthFit(i)=beta(5);
    widthErr(i)=err(1,5,1);
    %title(ha(i),sprintf('Wid %.0f ueV, lev %.2f',params(7)*1000,params(2)));    
end

tempFit=widthFit.*max(MCtemp*1000)./max(widthFit)
tempFitErr=widthErr.*max(MCtemp*1000)./max(widthFit)
figure(2); clf; errorbar(MCtemp*1000,tempFit,tempFitErr,'bo');
xlabel('MC temp (mK)'); ylabel('Temperature (mK)');
ax=axis;
axis([0 ax(2) 0 ax(4)]);
hold on;
[val ind]=max(widthFit);
plot([0 MCtemp(ind)]*1000,[0 tempFit(ind)],'r')
%%
%tempFit = widthFit.*max(MCtemp*1000)./max(widthFit)
figure(2); clf; 
%errorbar(MCtemp*1000,tempFit,tempFitErr,'bo');
plot(MCtemp,tempFit,'.-'); 
xlabel('MC temp (mK)'); ylabel('Temperature (mK)');
ax=axis;
axis([0 ax(2) 0 ax(4)]);
hold on;

[maxTemp,ind]=max(MCtemp); 
plot([0 maxTemp],[s0 tempFit(ind)],'r')