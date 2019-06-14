function atgradfix(bases,vstep,side)
% atgradfix(bases,dv,side)
%   Improves a specific set of basis vectors. By running charge 
%   scans, we see how those bases move the triple points, and add the X/Y
%   bases accordingly to stabilize. Iterate over increasingly large voltage changes,
%   as hopefully bases will move less and less. 
%   Later write in a way to average multiple runs. 
%
%   bases: can be 'all', 'leads', 'NTs', or any combination of gates. 
%   dv: the voltage step for each of the bases, e.g. atchg('N12',dv);
%   side: the tuneData set, or double dot to compute the basis for
%   Note: there are many functions to create gradients. For the x,y gates, use atxyfix2 (for chrg
%   scans), atxyfix3 (for stp/tl). For compensation matrix, use atcompfix3 (stp only). For other
%   gates, use atgradfix2 (chrg scans), atgradfix3 (stp/tl). 

global tuneData;

if ~exist('side','var') 
    sides ={tuneData.alternates.activeSetName};  
else 
    sides={side};
end
gateInds = 1:19; 
gtvals=cell2mat(smget(gateInds)); % MAKE ME SMARTER
[gates,b]=getgates(bases,sides); %puts selected bases into cell format. 

for i = 1:length(tuneData.alternates)
    autotune.swap(tuneData.alternates(i).activeSetName);
end

%For loop over the gates. For each gate, first run a charge scan to get original trip pts.
%Then increment by vstep, run scan for each side, and find how triple point
%moved. Update basis by adding in part of the X Y bases. 
for s=1:length(sides) %over sides    
    autotune.swap(sides{s});
    scan = smscanpar(tuneData.chrg.scan, tuneData.measPt);
    %scan.loops(1).rng = [-0.015 0.015]; scan.loops(2).rng = [-0.015 0.015];

    xy{1}=['X' upper(sides{s}(1))]; 
    xy{2}=['Y' upper(sides{s}(1))]; 
    d(1) = smrun(scan); 
    basxy=basislookup(xy);     
    for i=1:length(gates) 
        err=[];                
        tuneData.change(gates{i},vstep);
        fprintf('Changing %s by %.3f...\n',gates{i},vstep);
        d(i+1) = smrun(scan);
        tuneData.change(gates{i},-vstep)                
    end
end
smset(gateInds, gtvals);

for i = 1:length(d) 
    dataDiff = diff(d{i}, [], 2);        
    m=nanmean(dataDiff(:)); s=nanstd(dataDiff(:));
    dataDiff(abs(dataDiff)-m>3*s)=NaN;
    figure(70); clf;   
    imagesc(scan.loops(1).rng, scan.loops(2).rng, dataDiff);
    set(gca,'YDir','Normal')
    x=ginput(2);
    trip(i,1:2) = mean(x);
    if i > 1
        dtrippt = trip(i,:)-trip(1,:);
        grad = dtrippt'./vstep;             %Now we know how much things have moved.
        dbasis=tuneData.basis(1:4,basxy)*grad;        
        fprintf('Add X %.2f, Y %.2f \n to %s',grad(1), grad(2),gates{i-1});
        info=input('(y/n):', 's');
        doit = strcmp(info, 'y');
        if doit
            tuneData.basis(1:4,b(i-1))=tuneData.basis(1:4,b(i-1))-dbasis;
            tuneData.basisStore{tuneData.runNumber}=tuneData.basis; 
        end
    end
end
            
end

function [bases,b]=getgates(opts,sides)
bases=[];
if ~isempty(strfind(opts,'all')) 
    if length(sides)==2 
        bases={'Lead1', 'Lead2', 'Lead3', 'Lead4', 'T12', 'T34', 'N12', 'N34'};
    elseif length(sides)==1 && strcmp(sides,'right') 
        bases={'Lead3', 'Lead4', 'T34', 'N34','VRes'};
    elseif length(sides)==1 && strcmp(sides,'left') 
        bases={'Lead1', 'Lead2', 'T12', 'N12','VRes'};
    end
end
if ~isempty(strfind(opts,'leads'))
    bases={'Lead1', 'Lead2', 'Lead3', 'Lead4'};
end
if ~isempty(strfind(opts,'NTs')) 
    bases= {'T12', 'T34', 'N12', 'N34'};
end

if ~isempty(strfind(opts,'Lead1')), bases{end+1}='Lead1'; end  %#ok<*STREMP>
if ~isempty(strfind(opts,'Lead2')), bases{end+1}='Lead2';  end 
if ~isempty(strfind(opts,'Lead3')), bases{end+1}='Lead3';  end 
if ~isempty(strfind(opts,'Lead4')), bases{end+1}='Lead4';  end 
if ~isempty(strfind(opts,'T12')), bases{end+1}='T12';  end 
if ~isempty(strfind(opts,'T34')), bases{end+1}='T34';  end 
if ~isempty(strfind(opts,'N12')), bases{end+1}='N12';  end 
if ~isempty(strfind(opts,'N34')), bases{end+1}='N34';  end 
if ~isempty(strfind(opts,'VRes')), bases{end+1}='VRes';  end 

bases=unique(bases); 
b=basislookup(bases);
end                

function b=basislookup(basis)
global tuneData
%give this either a single gate in string form. 
if ischar(basis)
    basis={basis}; 
end
if iscell(basis)
    for j = 1:length(basis)
        b(j)=find(strncmp(basis{j},tuneData.baseNames,length(basis{j}))); 
    end
else 
    fprintf('Please input a cell or a string')
    b=nan; 
end
end