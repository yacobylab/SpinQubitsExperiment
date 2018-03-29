function [t1,ratio] =att1(side,scantime,opts,scan,data)
%[t1, ratio]= att1(side,scantime,opts,scan) Return t1 time appropriate to scantime.
%     side; left or right
% scantime; defaults to now
%     opts; 'before' t1 was measured before scan (default)
%           'after'  t1 was measured after scan.
%           'ask'    ask user.
if ~exist('opts','var'), opts='before'; end
if isopt(opts,'ask')
    [bt1,bt1r] = att1(side,scantime,'before');
    [at1,at1r] = att1(side,scantime,'after');
    if isnan(at1)
        t1=bt1;
        ratio=bt1r;
        return;
    end
    fprintf('T1 before: %g\n',bt1); 
    fprintf('T1 after: %g\n',at1);
    a=input('[B]efore or [a]fter?', 's');
    if ~isempty(a) && (a=='a' || a=='A')
        t1=at1;
        ratio=at1r;
        fprintf('Using after\n');
    else
        t1=bt1;
        ratio=bt1r;
        fprintf('Using before\n');
    end
    return
end
before = ~isopt(opts,'after');
if ~exist('scantime','var') && ~exist('scan','var')
    scantime = now; 
elseif ~exist('scantime','var') 
    scantime = getscantime(scan,data);
end
t1 = findt1(scantime,side,before);
if t1 < 1e-6
    warning('Short T1. Making a guess.');
    t1=20e-6;
end
if nargout > 1
    try
        %dict=scan.data.pulsegroups(1).dict{end};
        %mt=dict.meas(end); %end necessary for strange readout.
        %mt=mt.time(1)-(mt.time(4)+mt.time(5));
        measTime=scan.data.pulsegroups(1).readout.readout(3);
    catch
        warning('Guessing for measurement time.');
        measTime=2;
    end
    ratio=1e-6*measTime/t1;
end
end

function t1=findt1(scantime, side,before)
global tuneData;
t1=nan;
autotune.swap(side);
time = now; 
for i=length(tuneData.t1.t1):-1:1
    if ~isnan(tuneData.t1.t1(i))        
        time(end+1)=getFileTime(fullfile(tuneData.dir,sprintf('sm_t1%s_%04i.mat',upper(side(1)),i)));        
        if before && time(end)<scantime && time(end-1) > scantime 
            t1 = tuneData.t1.t1(i);
            break
        elseif ~before && time(end)<scantime && time(end-1) > scantime
            if ~exist('ind','var')
                t1 = tuneData.t1.t1(i);
            else
                t1 = tuneData.t1.t1(ind);
            end
            break
        end
        ind = i; 
    end
end
end