function flipJunc(opts)
% function will flip the direction of the junction scans.
% opts:
%   up: junction will flip to (2,0) being up and (1,1) down
%   down: junction will flip to (1,1) being up and (2,0) down
global tuneData
if ~exist('opts','var'), opts = ''; end
side = tuneData.activeSetName;
plsDict = pdload(side);

if isopt(opts,'up')
    tuneData.sepDir = [-1,1];
    tuneData.loadPos.slope = -2;
    tuneData.tl.slope = -2;
    plsDict.reload.val = [-2 -3];
    plsDict.rand(2).val = [-6 6];
    updateExch(struct('exch',[-1,1],'opts','all'))
end
if isopt(opts,'down')
    tuneData.sepDir = [1,-1];
    tuneData.loadPos.slope = -.5;
    tuneData.tl.slope = -.5;
    plsDict.reload.val = [-2 3];
    plsDict.rand(2).val = [6 -6];
    updateExch(struct('exch',[1,-1],'opts','all'))
end
pdsave(side,plsDict)
tuneData.updateAll('nodict');
fprintf('Flipping direction of junction... \n');
fprintf('Setting values, changing sep direction... \n');
fprintf('Random load position in (1,1) is now %d,%d \n',plsDict.rand(2).val(1),plsDict.rand(2).val(2));
fprintf('Load position is now %d,%d \n',plsDict.reload.val(1),plsDict.reload.val(2));
end