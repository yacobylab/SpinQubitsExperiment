function prettyFile = removePath(file)
lastDir = strfind(file,'\'); 
lastDir2 = strfind(file,'/'); 
if ~isempty(lastDir)
    file = file(lastDir(end)+1);
elseif ~isempty(lastDir2)
    file = file(lastDir2(end)+1);
end
prettyFile = file(1:end-4);
end