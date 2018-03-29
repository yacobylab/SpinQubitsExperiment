classdef (Abstract) Op < handle
    %Op < handle % generic autotune operation
    %   This class is really just a template for the different operations
    %   to add to autotune.Data
    %   inherited classes must overwrite the methods:
    %       getData: how to return data from a particular runNumber
    %       makeNewRun: initialize data when a new run is needed
    %       run: how to actually run a scan, take data, etc.
    
    %   Method overrideProperty will be inherited to set protected props
    
    properties (Abstract)
        subPlot;% must be assigned in each child class
    end
    
    %these are template methods. must be overridden
    %http://www.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html
    methods (Abstract) 
        getData(obj,runNumber);
        makeNewRun(obj,runNumber);
        run(obj,runNumber);
    end
    
    methods
        function overrideProperty(this,prop,val,opts)
            if ~exist('opts','var')
                opts = '';
            end
            if ~isprop(this,prop)
                fprintf('no property %s found',prop);
                return
            end
            q = sprintf('Override property %s with value with %g elements?. [y/n]  ',prop,numel(val));
            if ~isempty(strfind(opts,'noconfirm')) || strcmpi(input(q,'s'),'y' )
                this.(prop) = val;
            end
        end
    end
    
end

