function [slp,CBsep,msg]=anaCB
% Click on coulomb blockade plot to determine the separation and slope of
% lines.
% function [slp,CBsep,msg]=anaCB

sens = 1;
for i =1:2
    pkVals=ginput(2);
    yVals(i,:)=pkVals(:,2)';
    xVals(i,:)=pkVals(:,1)';
    if ~sens
        [~,xInd]=min(abs(xVals-xVals(i,1)));
        [~,yInd]=min(abs(yVals-yVals(i,1)));
        CBsz(i) = data(yInd,xInd);
    end
end

slp = diff(yVals,[],2)./diff(xVals,[],2);
b = yVals(:,1) - slp.*xVals(:,1);
CBsep = [diff(b), -diff(b)/slp(1)];

if ~sens
    Vout = 40e-6; I0 = 5.5e-9;
    R0 = Vout / I0;
    R = Vout / min(CBsz) - R0;
end
msg = sprintf('Slopes are is %3.3f and %3.3f y / x. Separation is %3.3f mV along Y, %3.3f mV along X  \n', slp(1),slp(2), 1e3*CBsep(1),1e3*CBsep(2));
fprintf('%s', msg);
end