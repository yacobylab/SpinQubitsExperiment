function printDict(side)
% Print current pulse dictionary on given side.
% function printDict(side)
% right, left
% In future, add  more documentation so that prints descriptive info about
% all elements. 
if ~exist('side','var')
    global tuneData; 
    side = tuneData.activeSetName;     
end
if ~isstruct(side)
    dict = pdload(side); 
else
    dict = side; 
end

fprintf('Single Element Dictionary Elements \n');
fprintf('%-11s %-11s %-22s %-20s \n','Name','Type','Time','Val')

dictElems = fieldnames(dict);
for i = 1:length(dictElems)
    currElem = dict.(dictElems{i});
    if numel(currElem) == 1 && ~isnumeric(currElem) % Check single element, has name. 
        times = sprintf('%3.3g ',currElem.time);
        eps = sprintf('%3.3g ',currElem.val);
        fprintf('%-11s %-11s %-22s %-20s \n',dictElems{i},currElem.type,times,eps);
    end
end
dictElemNames = {'wait','reload','adprep','meas'};
expl={'wait for time at val','ramp time, wait time at load, wait time at meas, load val','ramp from val 1 to 2 in time in dir val 3:4','wait for time at val'};
fprintf('\n')
fprintf('Multiple Element Dictionary Elements \n')
fprintf('%-9s %-10s %-10s \n','Type','Time','Val')
for i = 1:length(dictElems)
    currElem = dict.(dictElems{i});
    if numel(currElem) > 1 && ~isnumeric(currElem) % Check multi element, has name. 
        fprintf('%-11s: \n',dictElems{i})
        for j= 1:length(currElem)
            if isfield(currElem,'time')
                times = sprintf('%3.3g ',currElem(j).time);
            else
                times = '';
            end
            if isfield(currElem,'val')
                eps = sprintf('%3.3g ',currElem(j).val);
            else
                eps = '';
            end
            if ~isfield(currElem,'type'), currElem(j).type =''; end
            fprintf('%d %-11s %-10s %-10s \n',j,currElem(j).type,times,eps);
        end
    end
end
end