function makeNewTuneData(newDir)
clear global tuneData
global tuneData;
tuneData = autotune.Data;
tuneData.dir =[newDir '\tune_' datestr(now,'yyyy_mm_dd')];
% Create left as main, right as alternate. 
proplist  = {'chrg','zoom','line','stp','tl', 'tmp','loadTime','loadPos','lead','t1','twoSen'};
for j = 1:length(proplist)
   addprop(tuneData,proplist{j});
   tmp = proplist{j}; tmp(1) = upper(tmp(1));
   tuneData.(proplist{j}) = autotune.(tmp);
end
%tuneData.equalizeRuns; 

%tuneData.alternates = autotune.Data;
tuneData.alternates  = autotune.Data('right');
autotune.swap('right'); 
tuneData.addGroups; 
for j = 1:length(proplist)
    addprop(tuneData,proplist{j});
    tmp = proplist{j}; tmp(1) = upper(tmp(1));
    tuneData.(proplist{j}) = autotune.(tmp);
end
if ~isfolder(tuneData.dir), mkdir(tuneData.dir); end 
%% how to add a new field to tuneData. 
%newProp = 'twoSen';
%
%addprop(tuneData,newProp);
%tmp = newProp; tmp(1) = upper(tmp(1));
%tuneData.(newProp) = autotune.(tmp);
%addprop(tuneData.alternates(1),newProp);
%tuneData.alternates(1).(newProp)=autotune.(tmp);