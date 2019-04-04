function formatFig(n,opts,nrow,ncol)
% For all, rotate Y axes. Set axes and legend font sizes. 
% Quickly format figures to look nice without taking up too much space
% function formatFig(n,opts,nrow,ncol)
% n: A list of figure numbers to format.
% Sets font size, YTick rotation, position (i.e. size). 
% opts: 
% qpc2, qpc: for line scans. Set onscreen position differently, diff font size. 
% chrg, sens: color plots
% exch: 
% full: Compactify everything by moving labels closer to axis. 
% otherwise: 
% No option: set a line width. 

if ~exist('opts','var'), opts = ''; end
for i = 1:length(n)
    set(0,'CurrentFigure',n(i)); f = gcf; 
    axesInds = isgraphics(f.Children,'axes');
    legInds = isgraphics(f.Children,'legend');    
    if isopt(opts, 'qpc2')
        fontSize = 9; 
        tickRotation = -30; 
        legFontSize = 6; 
        f.Position = [746 9 932 977];
    elseif isopt(opts,'qpc') 
        fontSize = 7; 
        tickRotation = -30; 
        legFontSize = 5; 
        f.Position = [500 150 800 600];
    elseif isopt(opts,'exch')
        fontSize = 9; 
        tickRotation = -30; 
        legFontSize = 6; 
        f.Position = [500 150 800 600];
    elseif isopt(opts,'chrg')
        fontSize = 8.5; 
        tickRotation = -30; 
        legFontSize = 7; 
        f.Position = [100 150 800 630];
    elseif isopt(opts,'sens')
        fontSize = 8.5; 
        tickRotation = -30; 
        legFontSize = 7;                         
        f.Position = [900 150 800 630];    
    end
    [f.Children(axesInds).FontSize]=deal(fontSize);
    if any(legInds)            
            [f.Children(legInds).FontSize]=deal(legFontSize);
    end
    [f.Children(axesInds).YTickLabelRotation]=deal(tickRotation);
    if isopt(opts,'full') % This moves the X and Y labels closer to the axis. 
        axesNum = find(axesInds); 
        xLabOff = [32,16,10]; yLabOff = [23 11 7]; 
        for j = 1:length(axesNum)
            a = f.Children(axesNum(j)); 
            fitAxis(a);            
            a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/yLabOff(ncol);
            a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/xLabOff(nrow);            
        end       
    end
end
end