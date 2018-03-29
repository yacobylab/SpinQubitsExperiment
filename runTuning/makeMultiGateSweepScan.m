function scan=makeMultiGateSweepScan(scan,gateList,gateCount) 
% make scans that ramp different channels in different counts. 
% other channels are held at their initial value. 
% function scan=makeMultiGateScan(scan,gateList,gateCount) 
% send in 2d scan missing setchans, 2nd loop setchan / npoints 
% gateList: list of gates to scan 
% gateCount: cell with which point each channel is ramped in. 
scan.configfn(end+1).fn = @setTrafoArgs;
scan.configfn(end).args = {gateList}; 
maxGate=max(cell2mat(gateCount)); 
for i = 1:maxGate
    scan.data.setchan{i} = ''; 
end
for i = 1:length(gateList) 
    func = sprintf('@(x,y,p) p(1) + x(1)*(x(2)==%d',gateCount{i}(1)); 
    scan.data.setchan{gateCount{i}(1)} = [scan.data.setchan{gateCount{i}(1)} gateList{i} ', ']; 
    for j = 2:length(gateCount{i})
        func = [func sprintf('|x(2)==%d',gateCount{i}(j))];
        scan.data.setchan{gateCount{i}(j)} = [scan.data.setchan{gateCount{i}(j)} gateList{i}, ', ']; 
    end
    func = [func ')']; 
    scan.loops(1).trafofn(i).fn = str2func(func); 
end
scan.loops(2).setchan = 'count'; 
scan.loops(2).npoints = maxGate; 
scan.loops(1).setchan = gateList;
scan.disp(2) = scan.disp(1); 
scan.disp(2).dim = 2; 
for i = 1:maxGate 
    scan.data.setchan{i}(end-1:end)=[];
end
end