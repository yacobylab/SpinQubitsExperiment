function  val = awgcntrl(cntrl, chans)
% cntrl: stop, start, on off, wait, raw|amp, israw,  extoff|exton, isexton, err, clr, ison
% awgcntrl(cntrl, chans)
% start/stop : turn 'running' on and off. 
% on/off : turn voltage output of channels on and off. 
% wait: wait until all pending operations end (up to 600 s) 
% raw/amp: turn on raw (direct) mode, or amp (higher range, DAC channels added)
% exton/extoff: channels added from back on or off. 
% err: print out most recent (?) error. 
% clr: Clear all errors and print out number, last one. 
% norm / doubl: norm is range 0.6, double 1.2 . do we use this? 
% multiple commands given are processed in order.
% isamp and isexton return a vector of length chans specifying which are amp or in exton mode
% ison returns 0 and 1 for stopped and started, respectively, and .5 if waiting for trigger.

global awgdata;
val=[];
if ~exist('chans','var'), chans = []; end
breaks = [regexp(cntrl, '\<\w'); regexp(cntrl, '\w\>')]; % separate cntrl into separate words. 

for k = 1:size(breaks, 2)
    switch cntrl(breaks(1, k):breaks(2, k))
        case 'stop'
            for a=1:length(awgdata)
                fprintf(awgdata(a).awg, 'AWGC:STOP');
            end            
        case 'start'
            for a=1:length(awgdata)
                fprintf(awgdata(a).awg, 'AWGC:RUN');
            end
            awgcntrl('wait');            
        case 'off'
            for a=1:length(awgdata)
                for i = ch(awgdata(a), chans)
                    fprintf(awgdata(a).awg, 'OUTPUT%i:STAT 0', i);
                end
            end
        case 'on'
            for a=1:length(awgdata)
                for i = ch(awgdata(a), chans)
                    fprintf(awgdata(a).awg, 'OUTPUT%i:STAT 1', i);
                end
            end                        
        case 'wait'
            for a=1:length(awgdata)
                timeOut = awgdata(a).awg.timeout;
                awgdata(a).awg.timeout = 600;
                query(awgdata(a).awg, '*OPC?');
                awgdata(a).awg.timeout = timeOut;
            end
        case 'raw'
            if any(any(~awgcntrl('israw')))
                for a=1:length(awgdata)
                    if ~is7k(awgdata(a))
                        for i = ch(awgdata(a), chans)
                            fprintf(awgdata(a).awg, 'AWGC:DOUT%i:STAT 1', i);
                        end
                    end
                end
            else
                fprintf('Already raw\n');
            end            
        case 'amp'
            if any(any(awgcntrl('israw')))
                for a=1:length(awgdata)
                    if ~is7k(awgdata(a))
                        for i = ch(awgdata(a), chans)
                            fprintf(awgdata(a).awg, 'AWGC:DOUT%i:STAT 0', i);
                        end
                    end
                end
            else
                fprintf('Already amp\n');
            end            
        case 'israw'            
            val=[];
            for a=1:length(awgdata)
                if ~is7k(awgdata(a))
                    for i = ch(awgdata(a), chans)
                        fprintf(awgdata(a).awg, 'AWGC:DOUT%i:STAT?',i);
                        val(end+1) = fscanf(awgdata(a).awg,'%f');
                    end
                end
            end            
        case 'ison' %if instrument waiting for trigger, returns .5
            val=[];
            for a=1:length(awgdata)
                fprintf(awgdata(a).awg, 'AWGC:RST?');
                val(end+1) = .5*fscanf(awgdata(a).awg,'%f');
            end                        
        case 'exton'   %adds external DC to outputs specified in chans
            for a=1:length(awgdata)
                if ~is7k(awgdata(a))
                    for i = ch(awgdata(a),chans)
                        fprintf(awgdata(a).awg, 'SOUR%i:COMB:FEED "ESIG"', i);
                    end
                end
            end            
        case 'extoff'   %turns off external DC
            for a=1:length(awgdata)
                if ~is7k(awgdata(a))
                    for i = ch(awgdata(a),chans)
                        fprintf(awgdata(a).awg, 'SOUR%i:COMB:FEED ""', i);
                    end
                end
            end
        case 'isexton'
            val=[];
            for a=1:length(awgdata)
                if ~is7k(awgdata(a))
                    for i = ch(awgdata(a), chans)                        
                        fprintf(awgdata(a).awg, 'SOUR%i:COMB:FEED?\n',i);
                        val(end+1) = strcmp(fscanf(awgdata(a).awg, '%f'), 'ESIG');
                    end
                end
            end
        case 'err'
            for a=1:length(awgdata)
                err=query(awgdata(a).awg, 'SYST:ERR?');
                if ~strcmp(err(1:end-1), '0,"No error"'), fprintf('%d: %s\n',a,err); end
            end            
        case 'clr'
            for a=1:length(awgdata)
                i = 0;
                err2 = sprintf('n/a.\n');
                while 1
                    err = query(awgdata(a).awg, 'SYST:ERR?');
                    if strcmp(err(1:end-1), '0,"No error"')
                        if i > 0
                            fprintf('%d: %i errors. Last %s', a, i, err2);
                        end
                        break;
                    end
                    err2 = err;
                    i = i + 1;
                end
            end
        case 'norm'
            for i = 1:4
                fprintf(awgdata.awg, 'SOUR%i:VOLT:AMPL .6', i);
            end
        case 'dbl'
            for i = 1:4
                fprintf(awgdata.awg, 'SOUR%i:VOLT:AMPL 1.2', i);
            end
    end
end
end

function chans=ch(awg, chans)
if isempty(chans), chans=1:length(awg.chans); end
end

function val=is7k(awg)
val=length(awg.chans) <= 2;
end