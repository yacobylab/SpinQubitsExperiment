function [data,str] = formatData(data,cond,xvals,opts)
% function [data,str] = formatData(data,cond,xvals,opts)
% These are the functions called by options that format the data in
% plotChrgB. 
str='';
if isopt(opts,'glitch')
    m=nanmean(data(:)); s=nanstd(data(:));
    data(abs(data-m)>10*s)=NaN;
end
if isopt(opts,'cond')
    g0 = cond.i0 / cond.vout;
    gmeas = abs(data) / cond.vout;
    gsamp = (1 ./ gmeas - 1 / g0).^(-1);
    e= 1.6e-19; h = 6.6e-34;
    gqoc = 2 * e^2 / h;
    data = gsamp / gqoc;
end
if isopt(opts,'flat')
    [data,slp]=flattenData(data,xvals);
    for m = 1:size(slp,2)
        [n,x]=hist(slp(:,m));
        [~,histInd]=max(n);
        aveSlp(m) = x(histInd);
    end
    stdSlp = nanstd(slp);
    try             % We can use slp to characterize
        str = [str, sprintf('Sensor slopes at start, middle end: = %3.3f, %3.3f, %3.3f mV / mV. ',aveSlp(1),aveSlp(2),aveSlp(3))];
        str = [str, sprintf('Sensor std = %3.3f, %3.3f, %3.3f mV / mV \n',stdSlp(1),stdSlp(2),stdSlp(3))];
    end
    if isopt(opts,'glitch')
        m=nanmean(data(:)); s=nanstd(data(:));
        data(abs(data)-m>2*s)=NaN;
    end
end
if isopt(opts,'log')
    data=log10(data);
end
end