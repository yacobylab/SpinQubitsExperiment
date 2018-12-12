function replotRect(e,opts)
  if ~exist('opts','var'), opts = ''; end   
  ax = gca; 
  pl = ax.Children; 
  rect= getrect;
  xvals = pl.XData; yvals = pl.YData; 
  xInd = closePt(xvals,[rect(1),rect(1)+rect(3)]); 
  yInd = closePt(yvals,[rect(2),rect(2)+rect(4)]); 
  %r = rectangle('Parent',ax);
  %r.Position = rect;  
  if ~isopt(opts,'samefig') 
    figure(99); clf; 
  end
  data = pl.CData; 
  xInd = sort(xInd); yInd = sort(yInd); 
  imagesc(xvals(xInd(1):xInd(2)),yvals(yInd(1):yInd(2)),data(yInd(1):yInd(2),xInd(1):xInd(2))); 
  colorbar; 
  set(gca,'YDir','Normal')
  
end