function clim = setAllGraph(h)
% Return min and max color limits of a set of graphs. h is list of axes with color plots. 
% Won't crash if empty plots given, but may give meaningless data as empty
% plots and 2D plots have CLim [0,1]. 
% function clim = setAllGraph(h)
clim(1) = NaN;
clim(2) = NaN;
for i = 1:length(h)
    try
        clim(1) = min(clim(1), h(i).CLim(1));
        clim(2) = max(clim(2), h(i).CLim(2));
    catch
    end
end
end