function clearMask(chan)
global smdata; 

ic = smchaninst(chan); 
smdata.inst(ic(1)).cntrlfn([ic 6],[]);
end