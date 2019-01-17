function out=keyPressed(e,opts,file,j,cond)
% Function to add to figure as f.WindowKeyPressFcn, so that you can press
% key to run function on data in figure. 
%function out=keyPressed(e,opts,file,j,cond)
%   r: refresh file (useful when actively taking data) 
%	x: Load next set of scans in folder
%	s: Run spyview on data
%	m: measure distance between, slope of CB peaks
%	d: Click twice, gives distance between points.
%	p: Click once, give single point location
%	l: Drag a rectangle on screen and replot just that data in new figure. 
%   j: Drag a rectangle on screen and replot just that data in same figure.   
%   c: doublePt: Plot clicked point on both diff and charge (figures 1 and 3) 
%   e: hyst: change button down function to hyst.

switch e.Key
    case 'r'
        try
            plotChrg([opts ' last'],file,cond);
        catch
            warning('Couldn''t rerun program');
        end
    case 'm'
        [CBslp,CBsep,CBmsg]=anaCB;
        out.CBslp = CBslp;
        out.CBsep = CBsep;
        ppt = guidata(pptplot);
        msg = get(ppt.e_body,'String');
        CBlen = length(CBmsg);
        msgLen = size(msg,2);
        pad = msgLen- CBlen;
        padthing = '                                                                                                              ';
        if pad> 0
            CBmsg = [CBmsg, padthing(1:pad)];
        else
            CBmsg = [CBmsg(1:end+pad); CBmsg(end+pad+1:end), padthing(1:msgLen+pad)];
        end
        
        newmsg = [msg; CBmsg];
        set(ppt.e_body,'String',newmsg);
    case 's'
        smspyview(file{j})
    case 'x'
        plotChrg([opts, 'next'],file);
    case 'l'
        replotRect
    case 'j'
        replotRect([],'samefig')
    case 'd'
        d = ginput(2);
        xDist = d(2,1)-d(1,1);
        yDist = d(2,2)-d(1,2);
        fprintf('Distance is %3.3f mV in x dir and %3.3f mV in y dir \n',1e3*xDist,1e3*yDist);
    case 'p'
        x = ginput(1);
        fprintf('Point is at %3.3f V on x, %3.3f V on y \n',x(1),x(2));
    case 'c'
        f=e.Source;
        axesInds = find(isgraphics(f.Children,'axes'));
        for i =1:length(axesInds)
            f.Children(axesInds(i)).Children.ButtonDownFcn = @(src,clk) doublePt(src,clk);
        end
    case 'e'
        f=e.Source;
        axesInds = find(isgraphics(f.Children,'axes'));
        for i =1:length(axesInds)
            f.Children(axesInds(i)).Children.ButtonDownFcn = @(src,clk) hyst(src,clk);
        end
end
end