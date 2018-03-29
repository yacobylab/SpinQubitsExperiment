function scan = flipScan(scan)
% Used to flip the direction of the scan between each line. 
% function scan = flipScan(scan)

nsetchan1 = length(scan.loops(1).setchan); 
nsetchan2 = length(scan.loops(2).setchan);
scanCopy = scan;

scan.loops(1).setchan{1} = 'count'; 
scan.loops(1).trafofn(1).fn = @(x,y) x(1); % ????
scan.loops(1).trafofn(1).args = {}; 
for i = 1:nsetchan1
   scan.loops(1).setchan{i+1} = scanCopy.loops(1).setchan{i}; 
   if ~isfield(scanCopy.loops,'trafofn') || length(scanCopy.loops(1).trafofn)< i || isempty(scanCopy.loops(1).trafofn(i).fn)
        scan.loops(1).trafofn(i+1).args ={scanCopy.loops(1).rng,diff(scanCopy.loops(1).rng)./(scanCopy.loops(1).npoints-1)};
        scan.loops(1).trafofn(i+1).fn = @(x,y,rng,dr) rng(mod(x(2)+1,2)+1)-dr*(-1).^x(2)*(x(1)-1);     
   else
       scan.loops(1).trafofn(i+1).fn = @combineFunc; 
       scan.loops(1).trafofn(i+1).args{1} = @(x,y,rng,dr) rng(mod(x(2)+1,2)+1)-dr*(-1).^x(2)*(x(1)-1);     
       scan.loops(1).trafofn(i+1).args{2} = scanCopy.loops(1).trafofn(i).fn;
       scan.loops(1).trafofn(i+1).args{3} = {scanCopy.loops(1).rng,diff(scanCopy.loops(1).rng)./(scanCopy.loops(1).npoints-1)};
       scan.loops(1).trafofn(i+1).args{4} = scanCopy.loops(1).trafofn(i).args;
       scan.loops(1).trafofn(i+1).args{5} = 1; 
   end
end

scan.loops(1).rng = [1,scan.loops(1).npoints]; 
scan.loops(2).setchan{1} = 'count2';
scan.loops(2).trafofn(1).fn = @(x,y) x(2); 
scan.loops(2).trafofn(1).args = {};
for i = 1:nsetchan2
   scan.loops(2).setchan{i+1} = scanCopy.loops(2).setchan{i}; 
   if ~isfield(scanCopy.loops,'trafofn') || length(scanCopy.loops(2).trafofn)< i || isempty(scanCopy.loops(2).trafofn(i).fn)
        scan.loops(2).trafofn(i+1).args ={scanCopy.loops(2).rng,diff(scanCopy.loops(2).rng)./(scanCopy.loops(2).npoints-1)};
        scan.loops(2).trafofn(i+1).fn = @(x,y,rng2,dr2) rng2(1)+dr2*(x(2)-1);
   else
       scan.loops(2).trafofn(i+1).fn = @combineFunc; 
       scan.loops(2).trafofn(i+1).args{1} = @(x,y,rng2,dr2) rng2(1)+dr2*(x(2)-1);
       scan.loops(2).trafofn(i+1).args{2} = scanCopy.loops(2).trafofn(i).fn;
       scan.loops(2).trafofn(i+1).args{3} = {scanCopy.loops(2).rng,diff(scanCopy.loops(2).rng)./(scanCopy.loops(2).npoints-1)};
       scan.loops(2).trafofn(i+1).args{4} = scanCopy.loops(2).trafofn(i).args;
       scan.loops(2).trafofn(i+1).args{5} = 2; 
   end
end
scan.loops(2).rng = [1,scan.loops(2).npoints]; 
scan.data.flip = 1; 
end