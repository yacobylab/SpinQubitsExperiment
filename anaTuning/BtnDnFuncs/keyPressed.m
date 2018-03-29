function out=keyPressed(e,opts,file,j,cond)
%   r: refresh file
%	x: load next set of scans
%	s: spyview
%	m: measure distance between, slope of CB peaks
%	d: distance between two points.
%	p: single point location
%	l: replot rectange.
%   c: doublePt: plot clicked point on both diff and charge 
%   e: hyst: change button down function to hyst. 

switch e.Key
    case 'r'
        try
            plotChrgB([opts ' last'],file,cond);
        catch
            
            warning('Couldn''t rerun program');
        end
    case 'm'
        [CBslp,CBsep,CBmsg]=anaCB;
        out.CBslp = CBslp;
        out.CBsep = CBsep;
        ppt = guidata(pptplot3);
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
        plotChrgB([opts, 'next'],file);
    case 'l'
        replotRect(e)
    case 'd'
        d = ginput(2);
        xdist = d(2,1)-d(1,1);
        ydist = d(2,2)-d(1,2);
        fprintf('Distance is %3.3f mV in x dir and %3.3f mV in y dir \n',1e3*xdist,1e3*ydist);
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