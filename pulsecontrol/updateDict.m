function updateDict(opts,side) 
global tuneData; 
if ~exist('side','var') 
    side = tuneData.activeSetName;
end

if strcmpi(opts,'BR') 
    dir = -1; 
elseif strcmpi(opts,'TL')
    dir = 1; 
else 
    error('Invalid side. Use BR or TL'); 
end
dict = pdload(side); 
dict.exch.val = dir*[1,-1]; 
dict.sep.val = norm(dict.sep.val(1:2))*[1 -1]*dir/sqrt(2);
dict.rand(2).val(1:2) = dict.sep.val(1:2); 
sqrAmp = 4; 
if strcmpi(opts,'BR')
    leadSlp = [1 tuneData.chrg.yLeadSlope(end)]; % slope of vertical Left lead (lower left).
    leadSlp = leadSlp/norm(leadSlp);
    leadPos = -1e3*(tuneData.chrg.defaultOffset) + sqrAmp * leadSlp; %get the right offset
    tuneData.tl.slope = tuneData.chrg.yLeadSlope(end); 
else %y lead y square wave
    leadSlp = [1 tuneData.chrg.xLeadSlope(end)]; % slope of horizonal left lead. (upper left)
    leadSlp = leadSlp / norm(leadSlp);
    leadPos = -1e3*(tuneData.chrg.defaultOffset) - sqrAmp * leadSlp; %get the right offset
    tuneData.tl.slope = tuneData.chrg.xLeadSlope(end); 
end
dict.reload.val(1:2) = leadPos;

fprintf('Sep val set to %3.3f, %3.3f mV \n',dict.sep.val(1),dict.sep.val(2)); 
fprintf('Exch val set to %3.3f, %3.3f mV \n',dict.exch.val(1),dict.exch.val(2)); 
fprintf('Reload val set to %3.3f, %3.3f mV \n',dict.reload.val(1),dict.reload.val(2)); 
pdsave(side,dict);
end