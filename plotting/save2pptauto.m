function save2pptauto(text,figs)
% function save2pptauto(text,figs)
% text is a struct that contains all information for slide: 
% configch, configvals: saved channel values from smdata 
% consts: slides consts. 
% comments go under configvals. 
% if multiple images go on same slide, they will be tiled left to right,
% and slide width will increase to fit. all figure numbers in figs will be put on slide. 
% title 
% Adapated from save2word on Matlab File Exchange, created by
%Suresh E Joel, Mar 6,2003
%Virginia Commonwealth University
%Modification of 'saveppt' in Mathworks File Exchange
%and valuable suggestions by Mark W. Brown, mwbrown@ieee.org

global pptdata 
scale = 1; 
maxHeight = scale*525; maxWidth=scale*575;
wpicStart = scale*115; textWidth = scale*120; 
%% Setting up text 
if ~exist('text','var') || isempty(text) % Default title and body text:
  text = struct; 
end
text = def(text,'body',''); 
text = def(text,'consts',''); 
text = def(text,'opts','');
if isfield(text,'configch') && ~isempty(text.configch) % put configvals in struct text.consts.
    if isopt(text.opts, 'prevConfig')
        if ~isempty(pptdata.oldConfigvals)
            d.configch=text.configch; d.configvals=text.configvals; d.scan = text.scan; 
            chgStr=changeConfigGen(d,pptdata.oldConfigvals,pptdata.oldConfigch);
        end
        pptdata.oldConfigvals = text.configvals;
        pptdata.oldConfigch = text.configch;
        wpicStart = 10; % this may not always be true, may want consts + change
    else
        if isfield(text.scan,'consts') && ~isempty(text.scan.consts)
            vals = text.scan.consts;  %concatenate consts and vals
        else
            vals = struct;
        end
        for j =1:length(text.configch)
            vals(end+1).setchan = text.configch{j};
            vals(end).val = text.configvals(j);
        end
        text.consts=vals;
    end
else
    text.consts = '';
    wpicStart = 10; 
end
if isempty(text.consts) || isempty(fieldnames(text.consts)) % format configvals nicely, include time. 
    commentStr=''; %
else
    constStr=[];
    for i=1:length(text.consts) 
        if strcmp(text.consts(i).setchan,'Time') 
            constStr{i}=[text.consts(i).setchan,'     ',datestr(text.consts(i).val,'yy-mm-dd HH:MM:SS')];
        else
            constStr{i} = sprintf('%-12s %-4.4g', text.consts(i).setchan,text.consts(i).val);
        end
    end
    if ~any(strcmp(text.configch,'Time'))
        fileInd = find(strcmp({pptdata.dir.name},text.scanfile));
        if ~isempty(fileInd) 
            constStr{i+1} = ['Date: ' pptdata.dir(fileInd).date];
        end
    end
    commentStr=sprintf('%s\n',constStr{:});
end
if isfield(text,'comments') 
    commentStr=[commentStr text.comments];
end
if ~isempty(text.body) % convert text.body into correct formatting. 
    textCell=cellstr(text.body); % each row of text becomes cell. 
    bodyStr='';
    for i=1:length(textCell)
        bodyStr=[bodyStr textCell{i} newline];
    end
else
    bodyStr='';
end
if exist('chgStr','var'), bodyStr = [bodyStr,chgStr]; end
if isfield(text,'body2') && ~isempty(text.body2) % % convert text.body into correct formatting. 
    textCell=cellstr(text.body2);
    bodyStr2='';
    for i=1:length(textCell)
        bodyStr2=[bodyStr2 textCell{i} '\n'];
    end
end

%% PPT part
slideCount = get(pptdata.op.Slides,'Count'); % Get current number of slides
slideCount = int32(double(slideCount)+1); % Add a new slide (with title object):
newSlide = invoke(pptdata.op.Slides,'Add',slideCount,12); % Add new slide. 12 means blank slide. 
if ~isfield(text,'title') || isempty(text.title) % Insert text into the title object: FIXME
    picStart = scale*0; 
else
    set(newSlide.Shapes.Title.TextFrame.TextRange,'Text',text.title);    
    picStart = scale*80; 
end

pptdata.op.PageSetup.SlideHeight = scale*600; % Set slide height to ~ usual value
currWidth = pptdata.op.PageSetup.SlideWidth; 
textStart = 0;

pptdata.op.PageSetup.SlideWidth = max([currWidth, wpicStart+scale*85 + scale*550 * length(figs),1300]);  % Set width to be correct for max number of figures in single slide 
for i = 1:length(figs)    
    %set(0,'CurrentFigure',figs(i)); % FIX ME, not reallly working now so
    %changed to figure, which is slower. 
    figure(figs(i)); 
    print('-dmeta');  % copies current figure to clipboard
    %print('-dmeta'); 
    pic = invoke(newSlide.Shapes,'Paste'); % puts figure on new slide. 
    picHeight = get(pic,'Height'); picWidth = get(pic,'Width'); % Get height and width of picture:
    rat = picHeight / picWidth;
    if (picHeight/maxHeight > picWidth/maxWidth) % set height / width correctly to fit on slide. 
        set(pic,'Height',single(maxHeight));
        set(pic,'Width',single(maxHeight)/rat);
    else
        set(pic,'Width',single(maxWidth));
        set(pic,'Height',single(maxWidth)*rat);
    end
    set(pic,'Left',single(wpicStart)); % Center picture on right 3/4 of page (below title area):
    set(pic,'Top',single(picStart));
    wpicStart = wpicStart + get(pic,'Width'); % the next figure will be placed to right of this one. 
    textStart = max(textStart, get(pic,'Height')); % the text will go below the image. 
end
%% Text                                            % Left, Top, Width, Height
textConfig = invoke(newSlide.Shapes(1),'AddTextbox',1,scale*5,picStart,textWidth,scale*825); % to left, usually just configch. 
set(textConfig.TextFrame.TextRange.Font,'Size',scale*8);
set(textConfig.TextFrame.TextRange.Font,'Name','Consolas'); 
set(textConfig.TextFrame.TextRange,'Text',sprintf(commentStr));

textComments = invoke(newSlide.Shapes(1),'AddTextbox',1,textWidth + scale*30,textStart + picStart+8,scale*510,scale*175); % Below first image
set(textComments.TextFrame.TextRange.Font,'Size',scale*12);
set(textComments.TextFrame.TextRange,'Text',sprintf(bodyStr));

if exist('bodyStr2','var') % below second image. 
    textComments2 = invoke(newSlide.Shapes(1),'AddTextbox',1,textWidth + scale*30+scale*530,textStart + picStart+8,scale*650,scale*175);
    set(textComments2.TextFrame.TextRange.Font,'Size',scale*12);
    set(textComments2.TextFrame.TextRange,'Text',sprintf(bodyStr2));
end
end