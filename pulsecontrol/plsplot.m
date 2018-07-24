function plsplot(pulse, dict, ctrl)
% Plot pls using current dictionary elements. 
% plsplot(pulse, dict, ctrl)

if ~exist('ctrl','var')
    ctrl = '';
end
if ~exist('dict','var')
    dict = '';
end

styles={'b-','r-','g-'};
for i = 1:length(pulse)
    figure(30);
    if isempty(strfind(ctrl,'hold'))
        clf;
    else
        hold on;
        subplot(2,2,1); hold on;
        subplot(2,2,2); hold on;
        subplot(2,2,3); hold on;
    end
    pls = plstowf(pulse(i),dict);
    for c=1:length(pls.data)
        subplot(2,2,1); hold on;
        x=linspace(0,(size(pls.data(c).wf,2)-1)*1e9/pls.data(c).clk,size(pls.data(c).wf,2));
        plot(repmat(x,size(pls.data(c).wf,1),1)',pls.data(c).wf');        
        subplot(2,2,2); hold on;        
        for j = 1:2:size(pls.data(c).wf,1)-1
            plot(pls.data(c).wf(j, :), pls.data(c).wf(j+1, :),styles{(j+1)/2});
        end        
        if  any(pls.data(c).marker(:))
            subplot(2,2,3)
            plot(x,bitand(pls.data(c).marker, 1),x, bitand(pls.data(c).marker, 2)/2);
        end
    end
    if i < length(pulse) && isempty(strfind(ctrl, 'nopause'))
        pause;
    end
end
end