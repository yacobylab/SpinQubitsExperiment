function scan = sensorTrafofn( scan,fname,opts)
% Finds trafofn from sensor scan for sens scan. 
% function [ scan ] = sensorTrafofn( scan,fname,opts)
%   Opens a data file for a sensor scan, computes the appropriate trafofn, and modifies the input scan. 
%   You can 4 options: click once, and it will try to match that value with
%   line. Click twice, and it fits a line. Click > 2x and it will
%   interpolate between all points. If called with option 'auto', loads
%   most recent sensor scan and 
if ~exist('opts','var'), opts = ''; end
if isopt(opts,'auto')    
     fileList = dir; 
     fileNames = {fileList(:).name};
     sensList = strfind(fileNames,'sm_sensor');     
     sensInds = cellfun(@(x) ~isempty(x), sensList);
     sensNames = fileNames(sensInds); 
     dates = [fileList(sensInds).datenum]; % you may want to eliminate . and .. first.
     [~,inds] = sort(dates);
     file = sensNames(inds(end));      
     d = load(file{1}); 
else
    if ~exist('fname','var') || isempty(fname)
        file=get_files('sm_sensor*.mat');
        d=load(file{1});
    else
        d=load(fname);
    end
end
dataIndex=1;
data=d.data{dataIndex};
xvals = linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),d.scan.loops(1).npoints);
yvals = linspace(d.scan.loops(2).rng(1),d.scan.loops(2).rng(2),d.scan.loops(2).npoints);
xLab = d.scan.loops(1).setchan;
yLab = d.scan.loops(2).setchan;
if iscell(yLab) 
    yLab = char(yLab{:});
    yLab = yLab(1,:); 
end
if iscell(xLab) 
    xLab = char(xLab{:});
end

if length(scan.loops(2).setchan)~=2
	warning('Need to add a 2nd setchan')
end

f=figure(167); clf; f.Name = 'Sensor Trafofn';
if isopt(opts,'res') 
    div = scan.data.div;
    vout = scan.data.vout; 
    vout = vout * div; 
    r0 = 1e4; 
    data = vout ./ data - r0; 
end
dataDiff=diff(data,[],2);
dData = xvals(2)-xvals(1);
dataFilt = dataDiff; 
dataNorm = bsxfun(@rdivide,dataFilt,diff(xvals));    
xDiff = (xvals(1:end-1)+xvals(2:end))/2;
m=nanmean(dataNorm(:)); s=nanstd(dataNorm(:));

subplot(2,1,2); 
imagesc(xvals,yvals,data); colorbar; 
set(gca,'YDir','Normal'); hold on; 
xlabel(xLab); ylabel(yLab);

subplot(2,1,1)
imagesc(xDiff,yvals,dataNorm);  
set(gca,'YDir','Normal'); hold on; 
title('Charge sensing'); ylabel(yLab);
caxis([m-3*s,m+3*s]); colorbar; 
if isopt(opts,'auto')    
    scan = createPolyLine(dataNorm, yvals, xvals,xDiff,scan,dData,opts);    
else
    [X,Y]=ginput();
    if length(X) == 1 % try to maintain constant slope.
        [~,xInd] = min(abs(xDiff-X));
        [~,yInd] = min(abs(yvals-Y));
        setVal = dataNorm(yInd,xInd); % this is the value we're tryign to match.
        scan = createPolyLine(dataNorm, setVal, yvals, xvals,xDiff,scan,dData,opts);
    elseif length(X)==2 % fit a line between two points 
        slope=diff(Y)/diff(X);
        scan.loops(2).trafofn(2).fn=str2func('@(x,y,p)p(1)*(x(2)-p(2))+p(3)'); % p(3) gives value of sensor chan when a chan at p(2). p(1) gives slope of
        scan.loops(2).trafofn(2).args = {};
        scan.loops(2).trafofn(2).args={[slope, X(1), Y(1)]}; 
        fprintf('{[%.3f %.3f %.3f]}\n',slope,X(1),Y(1));
        figure(167);
        subplot(2,1,2); line(X,Y,'Color','r','LineWidth',2);
        subplot(2,1,1); line(X,Y,'Color','r','LineWidth',2);
    else % interpolate between n points. 
        xstring='[';
        ystring='[';
        for i=1:length(X)
            xstring=sprintf('%s %.3f ',xstring,X(i));
            ystring=sprintf('%s %.3f ',ystring,Y(i));
        end
        xstring=sprintf('%s]',xstring);
        ystring=sprintf('%s]',ystring);
        scan.loops(2).trafofn(2).fn=@(x,y,p,q) interp1(p,q,x(2),'linear','extrap'); % interp has args xin, yin , xextrap Here, xin is agate, yin is sensor gate.
        scan.loops(2).trafofn(2).args = {};
        scan.loops(2).trafofn(2).args{1}=X;
        scan.loops(2).trafofn(2).args{2}=Y;
        fprintf('fn=@(x,y,p,q) interp1(p,q,x(2))\nargs{1}=%s\nargs{2}=%s\n',xstring,ystring);
        rng = linspace(scan.loops(2).rng(1), scan.loops(2).rng(2), scan.loops(2).npoints);
        Yvals = interp1(scan.loops(2).trafofn(2).args{1},scan.loops(2).trafofn(2).args{2},rng,'linear','extrap');
        subplot(2,1,2); plot(rng,Yvals,'Color','r','LineWidth',2);
        subplot(2,1,1); plot(rng,Yvals,'Color','r','LineWidth',2);
    end
end
end

function scan = createPolyLine(dataNorm, yvals, xvals,xDiff,scan,dData,opts)
stepFunc = @(x1,x2,x) heaviside(x-x1)-heaviside(x-x2); % x > x1, x < x2 is 1.
if xvals(1)>xvals(end)
    funcT = @(x,p) (p(1)+p(2)*x)*heaviside(x-p(3)); 
    funcB = @(x,p) (p(1)+p(2)*x)*heaviside(p(4)-x);
else    
    funcT = @(x,p) (p(1)+p(2)*x)*heaviside(p(4)-x); 
    funcB = @(x,p) (p(1)+p(2)*x)*heaviside(x-p(3));
end
% or just find largest... 
[~,inds]=max(abs(dataNorm));%-setVal)); % find signal closest to setVal for each xval. 
yData = yvals(inds); % find sensor gate value for each.
for i = 1:length(xDiff) 
    SensSens(i) = dataNorm(inds(i),i); %sensor sensitivity 
end
figure(17); clf; hold on; 
plot(xDiff, SensSens); 
if length(xvals) > 200
    nLines = 12;
else
    nLines = 6;
end
inc = floor(length(xDiff)/nLines); % each line has inc points. 
for i = 1:nLines
    vals = ((i-1)*inc+1):i*inc; % set of inds for line. 
    if i==nLines, vals = ((i-1)*inc+1):length(xDiff); end 
    p{i}=robustfit(xDiff(vals),yData(vals)); % find slope / offset of each ilne. 
    endVals=sort(xDiff([vals(1),vals(end)])); % find start and end value for line. 
    args{i} = [p{i}(1) p{i}(2) endVals(1)-abs(dData) endVals(2)]; % args: slope, offset,more neg value, more pos value 
    if xvals(1) < xvals(2)&& i == nLines
        args{i}(end) = args{i}(end) + abs(dData);
    elseif xvals(1) > xvals(2) && i == 1
        args{i}(end) = args{i}(end) + abs(dData);
    end
end
func = @(x,p) (p(1)+p(2)*x)*stepFunc(p(3),p(4),x); % line with slp p1, offset p2, nonzero > p3, < p4. 

if length(xvals)>200 
    traf = @(x,y,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12) (funcT(x(2),p1)+func(x(2),p2)+func(x(2),p3)+func(x(2),p4)+func(x(2),p5)+func(x(2),p6)+func(x(2),p7)+func(x(2),p8)+func(x(2),p9)+func(x(2),p10)+func(x(2),p11)+funcB(x(2),p12));
else
    traf = @(x,y,p1,p2,p3,p4,p5,p6) (funcT(x(2),p1)+func(x(2),p2)+func(x(2),p3)+func(x(2),p4)+func(x(2),p5)+funcB(x(2),p6));
end
scan.loops(2).trafofn(2).fn = traf;
scan.loops(2).trafofn(2).args = args;
rng = linspace(scan.loops(2).rng(1), scan.loops(2).rng(2), scan.loops(2).npoints);
for i = 1:length(rng)
    yD(i) = traf([0 rng(i)],0,args{:});
end
figure(167);
subplot(2,1,2); plot(rng,yD,'Color','r','LineWidth',2);
subplot(2,1,1); plot(rng,yD,'Color','r','LineWidth',2);
if isopt(opts,'autox') 
     global scandata;
     slp = scandata.config.outSlp./scandata.config.inSlp;
     meanVal = mean(scandata.sens.loops(1).rng);
     sepX = @(x,p) p(1)*(x(1)-p(2));
     func = @(x,p,p0) (p(1)+p(2)*(x(2)+sepX(x,p0)))*stepFunc(p(3),p(4),x(2)); % line with slp p1, offset p2, nonzero > p3, < p4. 
     if xvals(1)>xvals(end)
         funcT = @(x,p,p0) (p(1)+p(2)*(x(2)+sepX(x,p0)))*heaviside(x(2)-p(3));
         funcB = @(x,p,p0) (p(1)+p(2)*(x(2)+sepX(x,p0)))*heaviside(p(4)-x(2));
     else
         funcT = @(x,p,p0) (p(1)+p(2)*(x(2)+sepX(x,p0)))*heaviside(p(4)-x(2));
         funcB = @(x,p,p0) (p(1)+p(2)*(x(2)+sepX(x,p0)))*heaviside(x(2)-p(3));
     end
     if length(xvals)>200
         traf = @(x,y,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12) (funcT(x,p1,p0)+func(x,p2,p0)+func(x,p3,p0)+func(x,p4,p0)+func(x,p5,p0)+func(x,p6,p0)+func(x,p7,p0)+func(x,p8,p0)+func(x,p9,p0)+func(x,p10,p0)+func(x,p11,p0)+funcB(x,p12,p0));
     else
         traf = @(x,y,p0,p1,p2,p3,p4,p5,p6) (funcT(x,p1,p0)+func(x,p2,p0)+func(x,p3,p0)+func(x,p4,p0)+func(x,p5,p0)+funcB(x,p6,p0));
     end
     args{1} = [slp meanVal];
     args(2:length(scan.loops(2).trafofn(2).args)+1) = scan.loops(2).trafofn(2).args;     
     scan.loops(1).trafofn(2).fn=traf; 
     scan.loops(1).trafofn(2).args = args;
     scan.loops(2).trafofn=[];
end
end