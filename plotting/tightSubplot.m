function ha=tightSubplot(vecSub,opts)
% Wrapper on tight_subplot to configure workable margins for different
% numbers of subplots.
% vecSub(1) = nrow , 2 = ncol
% Returns list of axes for subplots, counting left to right then top to
% bottom.
% Cutoffs are for 2, 4 columns/rows.
% The gap between plots and margins at edges are specified.
% Fixme: Check if axes are specfied correctly.
% opts: 
%   smart: 
%   nox: 
%   nolabely/x 
%   title
%   vert
if ~exist('opts','var'), opts = ''; end
if contains(opts,'smart')
    if vecSub<=2
        nrow =1; ncol = vecSub;
    elseif vecSub <= 4
        nrow = 2; ncol = 2;
    elseif vecSub<=6
        nrow = 2; ncol = 3;
    elseif vecSub <= 9
        nrow = 3; ncol = 3;
    elseif vecSub <= 12
        nrow = 3; ncol = 4; 
    elseif vecSub <= 16
        nrow = 4; ncol = 4; 
    else
        warning('Maximum elements currently 16. Fix me'); 
    end
elseif contains(opts,'vert')
    if vecSub <=4
        nrow = vecSub; ncol = 1;
    else
        ncol = 2; nrow = ceil(vecSub/ncol);
    end
else
    nrow = vecSub(1); ncol = vecSub(2);
end
if contains(opts,'nox')
    gapVert = 0.015;
else     
    if nrow <=2
        gapVert = 0.075;
    elseif nrow<=4
        gapVert = 0.065;
    else
        gapVert = 0.055;
    end
end
if isopt(opts,'title')
    gapVert = gapVert *1.35; 
end
if nrow > 1
    margVert = [0.1, 0.065]; % lower, upper
else
    margVert = [0.12, 0.065]; % lower, upper
end
if ncol > 1
    margHorz = [0.07, 0.06]; % left right
else
    margHorz = [0.09, 0.06]; % left right
end
if ncol <= 2
    gapHorz = 0.09;
elseif ncol<=4
    gapHorz = 0.065;
else
    gapHorz = 0.043;
end
if isopt(opts,'nolabely') 
    gapHorz = gapHorz * 0.75; 
end
if isopt(opts,'nolabelx') 
    gapVert = gapVert * 0.75; 
end
ha = tight_subplot(nrow, ncol, [gapVert,gapHorz], margVert, margHorz);
end