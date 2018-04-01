function ha=tightSubplot(vecSub)   
%vecSub(1) = nrow , 2 = ncol
nrow = vecSub(1); ncol = vecSub(2); 
    gapH = 0.063; 
    margH = [0.06 0.045]; 
    margW = [0.06 0.1]; 
     if nrow <= 2
         gapV = 0.12;   
     else
         gapV = 0.073; 
     end
ha = tight_subplot(nrow, ncol, [gapH gapV], margH, margW);
end