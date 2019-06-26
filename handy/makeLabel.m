function label = makeLabel(lab)
try
    if ischar(lab)
        lab = {lab};
    end
    if iscell(lab{1})
        lab = [lab{:}];
    end
    label = sprintf('%s, ',lab{:});
    label(end-1:end)=[];
catch
    label = '';
end
end
