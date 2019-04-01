function [t1,ratio] =att1(side,scantime,opts,scan,data)
% Return t1 time appropriate to scantime by looking through tuneData.
%[t1, ratio]= att1(side,scantime,opts,scan,data) 
%     side: left or right
% scantime: If not given but scan given, defaults to scan time. If neither
% given, defaults to now. 
%     opts: 'before' t1 was measured before scan (default)
%           'after'  t1 was measured after scan.
%           'ask'    Find both t1s and ask user which to use. 
% ratio: ratio of measurement time to t1 time. 
if ~exist('opts','var'), opts='before'; end
if isopt(opts,'ask')
    [bt1,bt1Rat] = att1(side,scantime,'before');
    [at1,at1Rat] = att1(side,scantime,'after');
    if isnan(at1)
        t1=bt1;
        ratio=bt1Rat;
        return;
    end
    fprintf('T1 before: %g\n',bt1); 
    fprintf('T1 after: %g\n',at1);
    a=input('[B]efore or [a]fter?', 's');
    if ~isempty(a) && (a=='a' || a=='A')
        t1=at1;
        ratio=at1Rat;
        fprintf('Using after\n');
    else
        t1=bt1;
        ratio=bt1Rat;
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
if t1 < 1e-6 % Under a microsecond suggests poor fit or bad data. 
    warning('Short T1. Making a guess.');
    t1=20e-6;
end
if nargout > 1
    if isfield(scan.data.pulsegroups(1),'readout') && isfield(scan.data.pulsegroups(1).readout,'readout') 
        measTime=scan.data.pulsegroups(1).readout.readout(3);
    else
        warning('Guessing for measurement time.');
        measTime=2;
    end
    ratio=1e-6*measTime/t1;
end
end

function t1=findt1(scantime, side,before)
% Crawl through tuneData to find t1 closest to scan time. 
global tuneData;
t1=nan;
if ~strcmp(tuneData.activeSetName,side)
    autotune.swap(side);
end
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