function run(cntl,varargin)
%function run(cntl,varargin): run stuff for autotune
%   cntl is a string of the operations to do
%   elements of cntl can be either properties of tuneData which have a run
%   method (i.e. inherit from autotune.Op) or can be methods of tuneData
%   itself. I.e. the function will either call tuneData.(op)(), or
%   tuneData.(op).run();
%   cntl can have lots of elements, which will be executed in order.
%   example: autotune.run('zoom') will call tuneData.zoom.run()
%   followed by tuneData.chrg.run();

global tuneData;

%allowed operations are probably the dynamic properties of tuneData 
opList = setdiff(properties(tuneData),properties('autotune.Data'));
opList = ([opList;methods(tuneData)]);
ops = (strsplit(cntl));
ops = intersect(ops,opList);
if isempty(ops)
   fprintf('no allowed operations found in control: %s\n',cntl);
   return
end

fprintf('tune run %i %s continues...',tuneData.runNumber,tuneData.activeSetName);
for op = 1:length(ops)
    if isa(tuneData.(ops{op}),'autotune.Op')
        tuneData.(ops{op}).run(tuneData.runNumber);
    elseif ismethod(tuneData,ops{op})
        tuneData.(ops{op})();
    else
       error('something has gone wrong. trying to autotune %s',ops{op}); 
    end
end

end