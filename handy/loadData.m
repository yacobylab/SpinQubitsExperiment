function d=loadData
[f,fpath] = uigetfile('','MultiSelect','on');
for i = 1:length(f) 
    d(i) = load(f{i}); 
end
end