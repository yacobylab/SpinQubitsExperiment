function good = isopt(opts, arg,varargin)
% check if something is an opt using strfind. 
%function good = isopt(opts, arg)
% good = ~isempty(strfind(opts,arg)); 
% if ~good
%     good = ~isempty(strfind(upper(opts),arg));
% end
% if ~good 
%     good = ~isempty(strfind(lower(opts),arg));
% end
if ~exist('arg','var') 
    error('No arg given'); 
end
optsList = strsplit(opts); 
    if any(strcmpi(optsList,arg))
        good = 1; 
    else
        good = 0; 
    end
end

