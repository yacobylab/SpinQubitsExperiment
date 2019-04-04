function fplot2ppt2(fig, scanfile, data)
% function val = fplot2ppt2(fig, scanfile, data)
% fig is figure number to be saved
% scanfile is datafile that generated figure
% data is struct with fields
% data.slidetitle
% data.comments: text to be added to slide under constants
% data.body: text to be added to slide under figure
% data.pptsavefile: ppt file to save it to, be careful of working directory
% opts can be 'noconfigvals', which will not include all of the congfigvals
% in the slide

if ~isempty(scanfile)
    load(scanfile, 'scan','configch','configvals'); %load all three separately, dont load data    
    data.configch = configch; 
    data.configvals = configvals; 
    data.scan = scan; 
else
    data.configch = [];
    data.scan = [];
end

data.scanfile=scanfile;

if exist('fig', 'var')
    save2pptman(data, fig);
else
    save2pptman(data);
end