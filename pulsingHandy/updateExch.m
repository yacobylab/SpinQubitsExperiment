function updateExch(config)
% function updateExch(config)
 % Load dictionary and update values, save. 
 % If config is string, converted to opts. 
 % config: 
 %    exch: Change exch to be [x(1),x(2)]*sqrt(2)/norm(x) The normalization
 %      is because magnitude usually for [1,-1]. 
 %    all: update the adprep, adread, sep directions to be those of exch. 
 %    sep: updates the norm of sep, but not the direction. 
 %    ad: update the adprep and adread final value to config.ad (length1). 
 %  opts: 
 %      perp: Use the most recent triple points to define exhange to be in perpendicular dir.  
global tuneData
dict = pdload(tuneData.activeSetName);  
if ischar(config), config= struct('opts',config); end
if ~isfield(config,'opts'), config.opts = ''; end
if isfield(config,'exch')
    newExch = config.exch/norm(config.exch)*sqrt(2); 
    dict.exch.val(1:2) = newExch; 
end
if isopt(config.opts,'perp')
    juncSlp = tuneData.chrg.trTriple(end,:) - tuneData.chrg.blTriple(end,:);
    juncSlp = juncSlp(2)/juncSlp(1); 
    epsSlp = [-1, 1./juncSlp];
    newExch = epsSlp/norm(epsSlp)*sqrt(2);
    dict.exch.val(1:2) = newExch;
    dict.sep.val(1:2) = dict.exch.val*norm(dict.sep.val)/sqrt(2);
    dict.adread.val(3:4) = dict.exch.val(1:2);
    dict.adprep.val(3:4) = dict.exch.val(1:2);
end
if isopt(config.opts,'all')
    dict.sep.val(1:2) = dict.exch.val*norm(dict.sep.val)/sqrt(2); 
    dict.adread.val(3:4) = dict.exch.val(1:2); 
    dict.adprep.val(3:4) = dict.exch.val(1:2);
end
if isfield(config,'sep')     
    dict.sep.val(1:2) = dict.sep.val*config.sep/norm(dict.sep.val)*sqrt(2); 
end
if isfield(config,'ad')     
    dict.adprep.val(2) = config.ad;
    dict.adread.val(2) = config.ad;    
end
pdsave(tuneData.activeSetName,dict); 
end