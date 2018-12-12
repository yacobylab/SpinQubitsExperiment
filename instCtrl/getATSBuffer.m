function [out,tm] = getATSBuffer(chan, samprate, npts, opts, downsamp)
% Acquire and output data from Alazar DAQ
% function out = getATSBuffer(chan, samprate, npts, opts)
% reads channel chan, at a sampling rate of samprate for npts
% opts can be 'hw' to use hardware trigger, default is software.
% Configures DAQ, then arms, then triggers, then reads.
% Assumes there is 1 DAQ on system and that it's name starts with ATS.
% Assumes hardware triggering done with first AWG5000.

global smdata;
if ~exist('opts','var'), opts = ''; end
if isnumeric(chan), chan = ['DAQ', num2str(chan)]; end
if ~ischar(chan), error('Chan must be number or channel name'); end
smset('samprate',samprate)
inst = sminstlookup('ATS');
channum = str2double(chan(end));

func = smdata.inst(inst).cntrlfn;
if ~exist('downsamp','var'), downsamp = 1; end
func([inst, channum, 5],npts,samprate/downsamp);
%
% if force
%   smdata.inst(inst).data.nrec = [0 0];
% end
func([inst, channum, 4]); % arm

if isopt(opts,'hw')    
    smatrigAWG('AWG1'); %trigger
else
    func([inst, channum, 3]); %software trigger
end
tic; 
out=smget(chan);
tm = toc; 
end

