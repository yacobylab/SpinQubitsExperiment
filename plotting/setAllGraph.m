function clim = setAllGraph(h) 
% set all the graphs to have the same color limits. 
    clim(1) = NaN; 
    clim(2) = NaN; 
    for i = 1:length(h) 
        try
        clim(1) = min(clim(1), h(i).CLim(1)); 
        clim(2) = max(clim(2), h(i).CLim(2)); 
        end
    end
end