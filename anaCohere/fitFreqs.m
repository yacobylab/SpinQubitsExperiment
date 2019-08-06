function [jFit,freq,data,out]=fitFreqs(files,opts)
if ~exist('opts','var'), opts = ''; end
if ~exist('files','var'), files = getFiles; end
a=procPlsData(files,'noplot');
for i = 1 : length(a)
    params = [];
    freq{i} = scanRng(a(i).scan,1);
    nQub = sum(contains(a(i).scan.loops(1).getchan,'DAQ'));
    
    for m = 1:nQub
        if all(size(a(i).data{1})>1)
            data{i,m} = squeeze(nanmean(a(i).data{m}));
        else
            data{i,m} = squeeze(a(i).data{m});
        end
        powInd = strcmp(a(i).configch,'RFpow3');
        if ~isempty(powInd)
            pow(i) = a(i).configvals(powInd);
        else
            pow(i) = NaN;
        end
        a(i).pow = pow(i);
        for j=1:length(freq{i})
            if all(isnan(data{i,m}(j,:)))
                jFit{i,m}(j:length(freq{i}))=nan;
                mse{i,m}(j:length(freq{i}))=nan;
                err(j:length(freq{i}),:)=nan;
                break
            end
            if isfield(a(i).scan.data,'FBset') && isnan(a(i).scan.data.FBset(j))
                jFit{i,m}(j) = nan;
            else
                [params{i,m}(j,:),~,~,~,mse{i}(j),err(j,:)]=fitosc(1:size(data{i,m},2),data{i,m}(j,:),'phase noplot',[6,Inf]);
                jFit{i,m}(j)=1e3*params{i,m}(j,4)./(2*pi);
            end
        end        
        jFit{i,m} = abs(jFit{i,m});
        jFit{i,m}(jFit{i,m}>1e4) = nan;
        jCorr{i,m} = jFit{i,m};
        jCorr{i,m}(jCorr{i,m}>5e2)=nan;
        jErr{i,m} = err(:,4)' * 1e3/(2 * pi);
        
        jCorr{i,m}(jErr{i,m}>1000)=nan;
        jErr{i,m}(jErr{i,m}>1000)=nan;
        jCorr{i,m} = removeOutliers(jFit{i,m});
        a(i).jErr{m} = jErr{i,m};
        a(i).jCorr{m} = jCorr{i,m};
    end
end
out.params = params; out.a = a; 
%%
figure(1202); clf;
extPlot=0;
if isopt(opts,'filter'), extPlot = extPlot + 1; end
if isopt(opts,'amp'), extPlot = extPlot+1; end
if isopt(opts,'t2s'), extPlot = extPlot + 1; end
if isopt(opts,'phase'), extPlot = extPlot + 1; end
ax = tightSubplot(2+extPlot,'vert nox');
for i =1:length(ax)
    hold(ax(i),'on');
    if i<2+extPlot, ax(i).XAxis.Visible='off'; end
end
figure(1203); clf;
ha = tightSubplot(length(a)*nQub,'smart');

inds = 1:length(a);
for i = inds
    for m = 1:nQub
        std = nanstd(jCorr{i,m});
        imagesc(ha(i),1:size(data{i,m},2),freq{i}*1e-9,data{i,m}); title(ha(i),num2str(pow(i)));
        ha(i).YDir = 'normal';
        ind = 1;
        if isopt(opts,'filter')
            ylabel(ax(ind),'J (MHz)');
            plot(ax(ind),freq{i}*1e-9,jCorr{i,m},'.-'); ind = ind+1;
        end
        if isopt(opts,'amp')
            ylabel(ax(ind),'Amplitude');
            plot(ax(ind),freq{i}*1e-9,params{i,m}(:,2),'.-');
            plot(ax(ind),freq{i}*1e-9,params{i,m}(:,1),'.-');
            ind = ind+1;
        end
        if isopt(opts,'phase')
            ylabel(ax(ind),'Phase');
            phs = mod(params{i,m}(:,3),pi);
            plot(ax(ind),freq{i}*1e-9,phs,'.-');
            ind = ind+1;
        end
        if isopt(opts,'t2s')
            ylabel(ax(ind),'T_2^*');
            plot(ax(ind),freq{i}*1e-9,1./params{i,m}(:,6),'.-'); ind = ind+1;
        end
        
        plot(ax(ind),freq{i}*1e-9,jFit{i,m},'.-','DisplayName',num2str(pow(i)));
        legInd = ind;
        ylabel(ax(ind),'J (MHz)'); ind = ind+1;
        
        plot(ax(ind),freq{i}*1e-9,jErr{i,m},'.-');
        ylabel(ax(ind),'J Error (MHz)');
    end
end
legend(ax(legInd),'show','location','best');
xlabel(ax(end),'Frequency (GHz)');
ppt = guidata(pptplot);
set(ppt.e_file,'String',a(1).filename);
set(ppt.e_figures,'String','[1202 1203]');
set(ppt.e_title,'String','Switch');
set(ppt.e_body,'String',char(files));
set(ppt.exported,'Value',0);
end