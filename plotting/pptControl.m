function pptControl(opts,filename,folder)
% opts: start, end, save, load
% will load the pptdata filename.
%function pptControl(opts,filename,folder)
% start: start a new ppt. 
% end: close ppt 
% save
% load
global pptdata
if exist('filename','var') && ~isempty(filename)
    pptdata.filename = filename;
end
if exist('folder','var') && ~isempty(folder) 
    pptdata.folder = folder; 
end
if isopt(opts,'start')
    pptdata.ppt = actxserver('PowerPoint.Application'); % Start an ActiveX session with PowerPoint:
    pptdata.op = invoke(pptdata.ppt.Presentations,'Add');   % Create new presentation:
elseif isopt(opts,'end')
    invoke(pptdata.op,'Close');
    %invoke(pptdata.ppt,'Quit');
elseif isopt(opts,'save')
    filespec = fullfile(pptdata.folder,[pptdata.filename,'.pptx']);
    if ~exist(filespec,'file') || pptdata.start
        invoke(pptdata.op,'SaveAs',filespec,11);   % Save file as new:
        pptdata.start =0; 
    else
        invoke(pptdata.op,'Save');   % Save existing file:
    end
elseif isopt(opts, 'load')
    filespec = fullfile(pptdata.folder,[pptdata.filename,'.pptx']);
    if ~exist(filespec,'file')
        pptControl('start')
        return
    end
    pptdata.ppt = actxserver('PowerPoint.Application'); % Start an ActiveX session with PowerPoint:
    numopen = pptdata.ppt.Presentations.count;
    isOpen = false;
    for i = 1:numopen
        if strcmp(filespec,pptdata.ppt.Presentations.Item(i).FullName)
            pptdata.op = pptdata.ppt.Presentations.Item(i);
            isOpen = true;
        end
    end
    if ~isOpen
        pptdata.op = pptdata.ppt.Presentations.Open(filespec,[],[],0);
    end
end
