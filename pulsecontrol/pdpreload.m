function [pd]=pdpreload(pd)
% replace any dictionaries in pd with the contents of the dictionary.  Judicious use can vastly speed up pulse dictionaries
% function [pd]=pdpreload(pd)

if iscell(pd)
   for i=1:length(pd)
     pd{i} = pdpreload(pd{i});
   end
   return;
end

if ischar(pd)
    pd=pdload(pd);
end
end