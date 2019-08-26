function addGroups(grps)
global fbdata; 
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
end