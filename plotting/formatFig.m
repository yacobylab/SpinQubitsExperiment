function formatFig(n,opts,nrow,ncol)
% Tries to make things look nice without taking up too much space. 
% function formatFig(n,opts,nrow,ncol)
% n a list of figure numbers to format.
% opts: 
% For all, rotate Y axes. Set axes and legend font sizes. 
% qpc2, qpc are line scans. Set onscreen position differently, diff font size. 
% chrg, sens: 
% exch
% full 
% otherwise: 
% No option: set a line width. 

if ~exist('opts','var'), opts = ''; end
for i = 1:length(n)
    set(0,'CurrentFigure',n(i)); f = gcf; 
    axesInds = isgraphics(f.Children,'axes');
    legInds = isgraphics(f.Children,'legend');    
    if isopt(opts, 'qpc2')
        [f.Children(axesInds).FontSize]=deal(9);
        if any(legInds)            
            [f.Children(legInds).FontSize]=deal(6);
        end
        [f.Children(axesInds).YTickLabelRotation]=deal(-30);
        f.Position = [746 9 932 977];
        %[f.Children(axesInds).FontSize]=deal(7);        
    elseif isopt(opts,'qpc') 
        [f.Children(axesInds).FontSize]=deal(7);
        if any(legInds)            
            [f.Children(legInds).FontSize]=deal(5);
        end
        [f.Children(axesInds).YTickLabelRotation]=deal(-30);
        f.Position = [500 150 800 600];
    elseif isopt(opts,'exch')
        [f.Children(axesInds).FontSize]=deal(9);
        if any(legInds)
            [f.Children(legInds).FontSize]=deal(6);
        end
        [f.Children(axesInds).YTickLabelRotation]=deal(-30);
        f.Position = [500 150 800 600];
    elseif isopt(opts,'chrg')
        [f.Children(axesInds).FontSize]=deal(8.5);
        if any(legInds)            
            [f.Children(legInds).FontSize]=deal(7);
        end
        [f.Children(axesInds).YTickLabelRotation]=deal(-30);
        f.Position = [100 150 800 630];
    elseif isopt(opts,'sens')
        [f.Children(axesInds).FontSize]=deal(8.5);
        if any(legInds)            
            [f.Children(legInds).FontSize]=deal(7);
        end
        [f.Children(axesInds).YTickLabelRotation]=deal(-30);
        f.Position = [900 150 800 630];
    else
        ax = gca;
        h = ax.Children;
        if ~isempty(h)
            [h.LineWidth] = deal(2.5);
        end
        ax.Box = 'off';
        ax.FontSize = 24;
        ax.FontName = 'Arial';
    end
    if isopt(opts,'full')
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