function chanDisp(channels,chanvals)
% make a figure 999 type image of the channel values when loading old data.
% function chanDisp(channels,chanvals)

nchans = length(channels); chanvals = num2cell(chanvals); 
screenData=get(0,'MonitorPositions');
f= figure(111); clf; f.MenuBar= 'none'; f.Name=  'Ana Channels';
% if size(screenData,1)>1 
%    nScr = 2; 
% else 
%     nScr = 1; 
% end
nScr = 1; 
f.Position = [10+screenData(nScr,1), screenData(nScr,4)-50-14*nchans, 220, 14*nchans+20];
%f.Position = [10, screenData(1,4)-50-14*nchans, 220, 14*nchans+20];
func = @(x,y) sprintf('%-20s %.5g', x,y); 
str = cellfun(func,channels,chanvals,'UniformOutput',false);
dispChan = uicontrol;     
dispChan.Style = 'text'; 
%dispChan.Position = [10+screenData(2,1), 10, 200, 14*nchans]; dispChan.HorizontalAlignment = 'Left';
%dispChan.Position = [10+screenData(2,1), 550, 200, 14*nchans]; dispChan.HorizontalAlignment = 'Left';
dispChan.Position = [10, 10, 200, 14*nchans]; dispChan.HorizontalAlignment = 'Left';
dispChan.String =  str;
dispChan.BackgroundColor= [.8 .8 .8];

