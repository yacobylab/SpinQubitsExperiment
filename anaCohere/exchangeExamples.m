global smdata; 
fileList = sortFiles(smdata.files.dir);
fpat = 'sm_RamseyL_(\d{4}).mat';
[fpatNames,fileNumsC]=regexp(fileList,fpat,'match','tokens');
%%
patFiles = [fpatNames{~cellfun(@isempty,fpatNames)}]; 
%%
fileList = sortFiles(smdata.files.dir);
fpat = 'sm_RamseyEL_(\d{4}).mat';
[fpatNames,fileNumsC]=regexp(fileList,fpat,'match','tokens');
patFiles = [fpatNames{~cellfun(@isempty,fpatNames)}]; 
%%
fitopts = 'residuals colorplot nocenter plotfit';
f = [11, 10, 401]; 
pptControl('start')
%%
for i = 124:length(patFiles)
    out = mfitEchoFish(patFiles{i},struct('mfitopts','none','grps',[4 Inf],'fignum',10,'opts',fitopts));
    if isempty(out) || ~isfield(out,'s'), continue; end
    str = sprintf('File %s. \n alpha = %3.3g. T_2 = %3.3g us. Amp = %3.3g. Fish Amp = %3.3g. T_2* = %3.3g ns', out.file(4:end-4),out.alpha, out.t2, out.amp,out.fishAmp,out.t2s);
    str2 = sprintf('T1 = %3.3g us. Peak separation %3.3g mV. Fidelity %3.3f%%.',out.s.t1*1e6,1e3*diff(out.s.meanvals),out.s.fidelity*100);
    slideInfo.body = str; 
    slideInfo.body2 = str2; 
    slideInfo.comments = ''; slideInfo.title = ''; scanfile=patFiles{i};
    slideInfo.scanfile = patFiles{i};
    save2pptauto(slideInfo,f)    
end
%%
fileList = sortFiles(smdata.files.dir);
fpat = 'sm_RamseyL_(\d{4}).mat';
[fpatNames,fileNumsC]=regexp(fileList,fpat,'match','tokens');
patFiles = [fpatNames{~cellfun(@isempty,fpatNames)}]; 
%%
pptControl('start')
for i = length(patFiles):-1:200
    anaExchange(patFiles{i},{'rng',[6 Inf],'spsize',[2,3],'opts', 'autoppt ramsey'});
end
%%
anaExchange(patFiles{i},{'rng',[3 Inf],'spsize',[2,3],'opts', 'autoppt ramsey'});