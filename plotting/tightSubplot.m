function ha=tightSubplot(vecSub)
% Wrapper on tight_subplot to configure workable margins for different
% numbers of subplots. 
% vecSub(1) = nrow , 2 = ncol 
% Returns list of axes for subplots, counting left to right then top to
% bottom. 
% Cutoffs are for 2, 4 columns/rows. 
% The gap between plots and margins at edges are specified. 
% FIxme: Check if axes are specfied correctly. 
nrow = vecSub(1); ncol = vecSub(2);
if nrow <=2
    gapVert = 0.075;
elseif nrow<=4
    gapVert = 0.065;
else
    gapVert = 0.055;
end
margVert = [0.08 0.065];
margHorz = [0.06 0.1];
if ncol <= 2
    gapHorz = 0.09;
elseif ncol<=4
    gapHorz = 0.055;
else
    gapHorz = 0.043;
end
ha = tight_subplot(nrow, ncol, [gapVert,gapHorz], margVert, margHorz);
end