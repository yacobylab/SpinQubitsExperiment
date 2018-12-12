function ha=tightSubplot(vecSub)   
%vecSub(1) = nrow , 2 = ncol
nrow = vecSub(1); ncol = vecSub(2); 
if ncol <=2     
    gapVert = 0.075;
elseif ncol<=4 
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