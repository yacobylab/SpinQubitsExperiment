function pulse = plstotab(pulse)
% Convert 'elem' to 'tab' pulse format.
% pulse = plstotab(pulse)
% In output pulse, the First row is time, and 2:3 are channel data, 4:5 IQ. 
% Integer pulse is taken from the database.
% Fields of pulse:
% name, xval, and taurc remain unchanged.
% format: must be 'elem' or 'tab'. Nothing is done for 'tab';
% data: struct array with fields type, time, val corresponding to pulse elements
%           to be concatenated.  type is a string specifying the type of element, 
%           val and time specicy pulse voltages and times, respectively.
%           Their meaning and the format of val depend on type.
% 
% Possible type strings and corresponding interpretation of val:
% raw: insert [time; val] into pulse table.
% mark: add time' to marktab
% fill: stretch this element to make the total pulse duration equal to time.
%       Idea for future development: allow several fills, each spreading a subset.
%       Would need a second element to flush previous fill, could be fill without time.
% wait: stay at val (row vector, one entry for each channel) for duration time.
%       If val has 3 entries, third is a scaling factor for the first two.
% reload: reload pulse at val (row vector, one entry for each channel).
%         time: [ramp time, wait time at load point, wait time at (0, 0) after load] 
% meas: measurement stage at [0, 0] for time(1), RF marker delayed by time(2) and off time(3) before end of the stage.
%       [time(2) is lead delay, time(3) is negative tail delay. 
%       Optional val(1) is the readout tag. If it is given and not nan, time 4 and 5 set its delays with the same convention as for the marker.
%       Optional val(2,3) moves the measurement point away from 0,0.  Makes meas_o obsolete.
% ramp: ramp to val (row vector, one entry for each channel) in time.  opt val(3) is multiplier
% comp: measurement compensation at val(1:2) (one for each channel) for duration time(1). 
%       Ramps voltage to target and back over time(2) and time(3) at the beginning and 
%       end of the stage, respectively. If length(val)>=4, val(3:4) are used as final value.
%       The compensation value could be determined automatically, but this feature is not 
%       implemented yet.
% adprep: adiabatic ramp along second diagonal (epsilon) from val(1) to
% val(2), ramp duration time. If val(3:4) given, they give directionof ramp
% (otherwise in sep dir). 
% adread: same, going the other way.
% RFburst: make AWG output sin wave: time is duration (in us).
%       val(1) = freq in MHz, val(2:3) = amplitude of x and y, 
%       optional val(4:5) = offset to 2 and 3
%       optional val(6:7) = phase for each channel
% markerburst: a wait with markers firing
%       time(1)=wait/marker time, val(1:3) where to wait, (just like exch)
%       time(2:3) = fire marker this time early, fire marker this long after exch
%       val(4:7)= which marker to fire, put in 0 or 1 for C1M1,C1M2,C2M1,C2M2
% IQburst: A wait with markers firing
%       takes 4 awg channels. first two are qubit channels, last two are iq
%       channels
%       time(1)=wait/marker time, val(1:3) where to wait for qubit, (just like exch)
%       val(4:5) amplitude of iq pulse
%       time(2:3) = fire marker this time early, fire marker this long after exch
%       val(6:9)= which marker to fire, put in 0 or 1 for C1M1,C1M2,C2M1,C2M2
% IQBurstSequence: A string of IQbursts
%       takes 4 awg channels. first two are qubit channels, last two are iq
%       channels
%       time(1:2)=fire marker this time early, fire marker this long after exch
%       time(3:end)=wait/marker times
%       val(1:2)=ch1 mult, ch2 mult
%       val(3)=eps
%       val(4:7)= which marker to fire, put in 0 or 1 for C1M1,C1M2,C2M1,C2M2
%       N=number of pulses
%       val(7+2n-1:7+2n)=IQ amp for gate n for n=1:N
%       val(end)=number of nonzero pulses. Used to make the for loop
%       faster
% IQBurstSequenceVarEps: A string of IQbursts
%       takes 4 awg channels. first two are qubit channels, last two are iq
%       channels
%       time(1:2)=fire marker this time early, fire marker this long after exch
%       time(3:end)=wait/marker times
%       val(1:2)=ch1 mult, ch2 mult
%       val(3:6)= which marker to fire, put in 0 or 1 for C1M1,C1M2,C2M1,C2M2
%       N=number of pulses
%       val(6+3n-2,6+3n-1,6+3n)=IQ amps and eps for gate n for n=1:N
%       val(end)=number of nonzero pulses. Used to make the for loop
%       faster
% XYburst: A set of pi pulses alternately around Z and X axes. 
%       time(1) = length of each jpi pulse.
%       time(2) = length of each dbz pi pulse.
%       time(3) = total length of pulse. If not a multiple of time(1), new pulse added with < pi time. 
%       val(1) val of eps for dBz osc. 
%       val(2) multiplier for second channel, i.e eps_chan2 = eps_chan1 *
%       mult
%       val(3:end) vector of eps on chan1 for J pi pulses. Could be changed to function
%       later. 
% XYburst2: allows different times for dbz and jpi pulses
% XYburst3: allows different time and different epsilon values for dbz and jpi
% MARKTAB:
%  each column gives a marker pulse, length of time it is on. 
%  [ start_time; ch1_mk1_wid ; ch1_mk2_wid ; ch_2_mk1_wid ...]
% ie,
%  [ 2 ; 0 ; 1 ; 0 ; 1 ] fires markers ch1mk1 and ch2mk2 for 1 us starting at 2us.

pulse = plsdefault(pulse);
dt=-1e-9;  % Shortest meaningful length
warning('off','MATLAB:catenate:DimensionMismatch');
switch pulse.format
    case 'tab'
        return;      
    case 'elem'        
        pulsetab = zeros(5, 0);
        marktab =  zeros(9, 0);        
        comppos = []; fillpos = []; readout = []; readpos = [];        
        pulseDef = pulse.data;        
        for i = 1:length(pulseDef)
            switch pulseDef(i).type
                case 'raw'
                    pulsetab = [pulsetab, [pulseDef(i).time; pulseDef(i).val]];
                case 'mark'
                    marktab = [marktab, pulseDef(i).time'];
                case 'fill'
                    fillpos = size(pulsetab, 2);
                    filltime = pulseDef(i).time(1);
                    fillmarkpos = size(marktab,2);                    
                case 'wait'
                    if pulseDef(i).time(1) > 1e-11
                        fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        if length(pulseDef(i).val) > 2
                          pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(3)*pulseDef(i).val(1:2)', 1, 2);
                        else
                          pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(1:2)', 1, 2);
                        end
                    end
                case 'reload'
                    % If we're filling the load, push the fillpos 1 forward
                    % so we stretch the wait at the loadpos, not the ramp to the loadpos                    
                    % Ignore zero length loads
                    if pulseDef(i).time(2) > 1e-11
                      fillload = (fillpos == size(pulsetab,2));                        
                      pulsetab(1, end+(1:4)) = pulsetab(1, end) + cumsum(pulseDef(i).time([1 2 1 3]));
                      pulsetab(2:3, end+(-3:0)) = [repmat(pulseDef(i).val(1:2)', 1, 2), zeros(2)];
                      fillpos = fillpos + fillload;                    
                    end               
                case 'meas'
                    if length(pulseDef(i).val) == 3
                        measPt = pulseDef(i).val(2:3);
                    else
                        measPt = [0,0];
                    end
                    measStart = pulsetab(1, end) + dt;
                    measEnd = pulsetab(1, end) + pulseDef(i).time(1);
                    pulsetab(:, end+1) = [measStart; measPt'];   
                    pulsetab(:, end+1) = [measEnd; measPt'];                     
                    marktab(1,end+1) = pulsetab(1, end-2)+pulseDef(i).time(2);
                    markLength = pulseDef(i).time(1:3)*[1; -1; -1];
                    marktab(2:5, end) = [0; 0; 0; markLength];
                    if ~isempty(pulseDef(i).val) && ~isnan(pulseDef(i).val(1))                    
                        measTime = pulseDef(i).time([1 4 5])*[1; -1; -1];
                        useDelay = pulseDef(i).val(1);
                        %               flag for a delay on readout, the measurement start time and measurement duration    
                        readout(end+1, :) = [useDelay, measStart + pulseDef(i).time(4), measTime];
                        readpos(end+1) = size(pulsetab, 2)-2; % Which index readout starts at. 
                    end
                case 'ramp'  %allow for multiplies in ramps - helps get direction right
                    if length(pulseDef(i).val) ==3
                        mult = pulseDef(i).val(3);
                    else
                        mult = 1;
                    end
                    pulsetab(1, end+1) = pulsetab(1, end) + pulseDef(i).time(1);
                    pulsetab(2:3, end) = mult*pulseDef(i).val(1:2);
                case 'comp'
                    comppos = size(pulsetab, 2)+1;
                    compval  = pulseDef(i).val(1:2);
                    pulsetab(1, end+(1:4)) = pulsetab(1, end) + [0 pulseDef(i).time(2), pulseDef(i).time(1)-sum(pulseDef(i).time(2:3)), ...
                        pulseDef(i).time(1)];
                    pulsetab(2:3, end+(-3:0)) = 0;
                    if length(pulseDef(i).val) >= 4
                        pulsetab(2:3, end) = pulseDef(i).val(3:4);
                    end
                case 'adprep'
                    if pulseDef(i).time(1) > 1e-11
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(1)];
                        if(length(pulseDef(i).val) <= 2)
                            dir=[-1 1];
                        else
                            dir = pulseDef(i).val(3:4);
                        end
                        pulsetab(2:3, end-1) = pulseDef(i).val(1)  * dir;
                        pulsetab(2:3, end) = pulseDef(i).val(2) * dir;
                    end
                case 'adread'
                    if pulseDef(i).time(1) > 1e-11
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(1)];
                        if(length(pulseDef(i).val) <= 2)
                            dir=[-1 1];
                        else
                            dir = pulseDef(i).val(3:4);
                        end
                        pulsetab(2:3, end-1) = pulseDef(i).val(2)  * dir;
                        pulsetab(2:3, end) = pulseDef(i).val(1)  * dir;
                    end
                case 'RFburst' 
                    a = pulseDef(i).val;
                    a(isnan(a))=0;
                    a(1+length(a):7)=0;
                    offst = a(4:5);
                    ph = a(6:7);
                    ts = 0:abs(1e6*dt):pulseDef(i).time(1);
                    if ~exist('firstburst','var') % the first burst will determine the phase for the rest of the pulse
                       firstburst = pulsetab(1,end); 
                    end
                    if length(ts) > 1 % number of time points
                        if pulseDef(i).val(1) > 1/abs(dt) % more than 1GHz
                            warning('asking for frequency above 1GHZ');
                        end
                        t_off = pulsetab(1,end)-firstburst;                        
                        pulsetab(1,end+(1:length(ts)))=pulsetab(1,end)+ts;
                        pulsetab(2,end+(-length(ts)+1:0)) = offst(1)+pulseDef(i).val(2)*sin(ph(1)+2*pi*pulseDef(i).val(1)*(ts+t_off));
                        pulsetab(3,end+(-length(ts)+1:0)) = offst(2)+pulseDef(i).val(3)*sin(ph(2)+2*pi*pulseDef(i).val(1)*(ts+t_off));
                    else
                        warning('asking for only one point of oscillations');
                    end
                case 'markerburst'
                    ctime = pulsetab(1,end); %current time                    
                    if pulseDef(i).time(1) > 1e-11 %first treat it as a wait
                        fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        if length(pulseDef(i).val) > 2
                            pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(3)*pulseDef(i).val(1:2)', 1, 2);
                        else
                            pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(1:2)', 1, 2);
                        end                        
                        if length(pulseDef(i).time)>1 && all(~isnan(pulseDef(i).time(2:3))) % now deal with the markers
                           stdly = pulseDef(i).time(2); %start delay
                           edly = pulseDef(i).time(3); %end delay
                        else
                            stdly = 0; edly = 0;
                        end
                        mk = pulseDef(i).val(4:7); mk(isnan(mk))=0;
                        marktab = [marktab, [ctime-stdly, (pulseDef(i).time(1)+edly+stdly)*mk]'];
                    end                    
                case 'IQburst'
                    ctime = pulsetab(1,end); %current time
                    %first treat it as a wait
                    if pulseDef(i).time(1) > 1e-11
                        fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        
                        %Qubit
                        pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(3)*pulseDef(i).val(1:2)', 1, 2);

                        %IQ
                        pulsetab(4:5, end+(-1:0)) = repmat(pulseDef(i).val(4:5)', 1, 2);

                        % now deal with the markers
                        %marktab = [marktab, pulsedef(i).time'];
                        if length(pulseDef(i).time)>1 && all(~isnan(pulseDef(i).time(2:3)))
                            stdly = pulseDef(i).time(2); %start delay
                            edly = pulseDef(i).time(3); %end delay
                        else
                            stdly = 0; edly = 0;
                        end
                        mk = pulseDef(i).val(6:9); mk(isnan(mk))=0;
                        marktab = [marktab, [ctime-stdly, 0,0,0,0,(pulseDef(i).time(1)+edly+stdly)*mk]'];
                    end                                    
                case 'IQBurstSequence'
                    stdly = pulseDef(i).time(1); %start delay
                    edly = pulseDef(i).time(2); %end delay
                    for o=1:pulseDef(i).val(end)%for all the real pulses
                        ctime = pulsetab(1,end); %current time
                        %first treat it as a wait
                        if pulseDef(i).time(2+o) > 1e-11
                            fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                            %I dont even think fillpos is used in iqburst
                            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(o+2)]; %pinf.tbase*1e6/pinf.clk.                            
                            %Qubit
                            pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(3)*pulseDef(i).val(1:2)', 1, 2);%eps is fixed to be the same for all pulses                            
                            %IQ
                            pulsetab(4:5, end+(-1:0)) = repmat(pulseDef(i).val((7+2*o-1):(7+2*o))', 1, 2);%this is the huge list of iq pulse values                            
                            % now deal with the markers
                            mk = pulseDef(i).val(4:7); mk(isnan(mk))=0;
                            if isempty(marktab)
                                marktab = [marktab, [ctime-stdly, 0,0,0,0,(pulseDef(i).time(2+o)+edly+stdly)*mk]'];
                                stTime=ctime-stdly;
                            else
                                marktab(:,end) = [stTime, 0,0,0,0,(ctime+pulseDef(i).time(2+o)+edly-stTime)*mk]';
                            end
                        end
                    end                    
                case 'IQBurstSequenceVarEps'
                    stdly = pulseDef(i).time(1); %start delay
                    edly = pulseDef(i).time(2); %end delay
                    for o=1:pulseDef(i).val(end)%for all the real pulses
                        ctime = pulsetab(1,end); %current time
                        %first treat it as a wait
                        if pulseDef(i).time(2+o) > 1e-11
                            fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                            %I dont even think fillpos is used in iqburst
                            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulseDef(i).time(o+2)]; %pinf.tbase*1e6/pinf.clk.
                            
                            %Qubit
                            pulsetab(2:3, end+(-1:0)) = repmat(pulseDef(i).val(6+3*o)*pulseDef(i).val(1:2)', 1, 2);%eps is independent for all pulses                            
                            %IQ
                            pulsetab(4:5, end+(-1:0)) = repmat(pulseDef(i).val((6+3*o-2):(6+3*o-1))', 1, 2);%this is the huge list of iq pulse values                            
                            % now deal with the markers
                            mk = pulseDef(i).val(3:6); mk(isnan(mk))=0;
% %                             If things go bad, this next line is the one to uncomment
%                             marktab = [marktab, [ctime-stdly, 0,0,0,0,(pulsedef(i).time(2+o)+edly+stdly)*mk]'];
                            if isempty(marktab)
                                marktab = [marktab, [ctime-stdly, 0,0,0,0,(pulseDef(i).time(2+o)+edly+stdly)*mk]'];
                                stTime=ctime-stdly;
                            else
                                marktab(:,end) = [stTime, 0,0,0,0,(ctime+pulseDef(i).time(2+o)+edly-stTime)*mk]';
                            end
                        end
                    end 
                case 'XYburst' % factors needed below to avoid the pathological situation where matlab thinks that ceil(.07/.01)=8, but ceil(70/10)=7;
                    npulse = ceil((pulseDef(i).time(2).*1000)/(pulseDef(i).time(1).*1000));                    
                    finalPulseTime = mod(pulseDef(i).time(2),pulseDef(i).time(1));                     
                    mult = pulseDef(i).val(2); 
                    jEpsVals = pulseDef(i).val(3:end);                     
                    if length(jEpsVals) < ceil(npulse/2)
                        jEpsVals(end+1:ceil(npulse/2)) = jEpsVals(end); 
                    end
                    for k = 1:npulse
                        pulsetab(1,end+1) = pulsetab(1,end);%-dt;
                        pulsetab(1,end+1) = pulsetab(1,end-1) + pulseDef(i).time(1); 
                        if mod(k,2)                            
                            pulsetab(2:3,end-1:end) = repmat(jEpsVals((k+1)/2).*[1; mult],[1,2]);
                        else
                            pulsetab(2:3,end-1:end) = repmat(pulseDef(i).val(1) .* [1; mult],[1,2]); 
                        end                        
                    end
                    if finalPulseTime ~= 0 
                        pulsetab(1,end) = pulsetab(1,end-2)+finalPulseTime; 
                    end                    
                case 'XYburst2'
                    %allows different times for dbz and J rotations                    
                    %npulse now means how many combined (jpi dbz) pulses there are.
                    npulse = floor((pulseDef(i).time(3).*1000)/(sum(pulseDef(i).time(1:2)).*1000));
                    njpi=npulse;
                    ndbz=npulse;
                    
                    finalPulseTime = mod(pulseDef(i).time(3),sum(pulseDef(i).time(1:2)));
                    if finalPulseTime>0 && finalPulseTime<=pulseDef(i).time(1)
                        njpi=njpi+1;
                    elseif finalPulseTime>0
                        njpi=njpi+1;
                        ndbz=ndbz+1;
                        finalPulseTime=finalPulseTime-pulseDef(i).time(1);
                    end                    
                    mult = pulseDef(i).val(2);
                    jEpsVals = pulseDef(i).val(3:end);
                    if length(jEpsVals) < njpi
                        jEpsVals(end+1:njpi) = jEpsVals(end);
                    end
                    for k = 1:njpi+ndbz
                        pulsetab(1,end+1) = pulsetab(1,end);%-dt;
                        pulsetab(1,end+1) = pulsetab(1,end-1) + pulseDef(i).time(2-mod(k,2)); %switch times depending if jpi or dbzpi
                        if mod(k,2)
                            pulsetab(2:3,end-1:end) = repmat(jEpsVals((k+1)/2).*[1; mult],[1,2]);
                        else
                            pulsetab(2:3,end-1:end) = repmat(pulseDef(i).val(1) .* [1; mult],[1,2]);
                        end
                    end
                    if finalPulseTime ~= 0
                        pulsetab(1,end) = pulsetab(1,end-1)+finalPulseTime;
                    end                    
                case 'XYburst3'%allows different times for dbz and J rotations, different eps values for dbz.
                    %npulse now means how many combined (jpi dbz) pulses there are.
                    npulse = floor(round(pulseDef(i).time(3).*1e7)/round(sum(pulseDef(i).time(1:2)).*1e7));
                    njpi=npulse;
                    ndbz=npulse;                    
                    finalPulseTime = mod(pulseDef(i).time(3),sum(pulseDef(i).time(1:2)));
                    if finalPulseTime>0 && finalPulseTime<=pulseDef(i).time(1)
                        njpi=njpi+1;
                    elseif finalPulseTime>0
                        njpi=njpi+1;
                        ndbz=ndbz+1;
                        finalPulseTime=finalPulseTime-pulseDef(i).time(1);
                    end                    
                    mult = pulseDef(i).val(1);
                    epsVals = pulseDef(i).val(2:end);
                    neps=length(epsVals);
                    jEpsVals=epsVals(1:2:neps);
                    dbzEpsVals=epsVals(2:2:neps);                    
                    if length(jEpsVals) < njpi
                        jEpsVals(end+1:njpi) = jEpsVals(end);
                    end                    
                    if length(dbzEpsVals) < ndbz
                        dbzEpsVals(end+1:ndbz) = dbzEpsVals(end);
                    end                    
                    for k = 1:njpi+ndbz
                        pulsetab(1,end+1) = pulsetab(1,end);%-dt;
                        pulsetab(1,end+1) = pulsetab(1,end-1) + pulseDef(i).time(2-mod(k,2)); %switch times depending if jpi or dbzpi
                        if mod(k,2)
                            pulsetab(2:3,end-1:end) = repmat(jEpsVals((k+1)/2).*[1; mult],[1,2]);
                        else
                            pulsetab(2:3,end-1:end) = repmat(dbzEpsVals((k)/2).*[1; mult],[1,2]);
                        end
                    end
                    if finalPulseTime ~= 0
                        pulsetab(1,end) = pulsetab(1,end-1)+finalPulseTime;
                    end                    
                otherwise
                    error('Invalid pulse element %i: %s.\n', i, pulseDef(i).type)
            end
        end
        if ~isempty(comppos)
            pulsetab(2:3, comppos+(1:2)) = 2*repmat(compval(1:2)', 1, 2);
        end
        if ~isempty(fillpos)
            filltime = filltime - pulsetab(1, end);
            if filltime < 0
                error('Pulse too long by %g (target %g).',-filltime,filltime+pulsetab(1,end));
            end
            pulsetab(1, fillpos+1:end) = pulsetab(1, fillpos+1:end) + filltime;
            if ~isempty(readpos)
                readout(readpos > fillpos, 2) = readout(readpos > fillpos, 2) + filltime;
            end
            marktab(1, fillmarkpos+1:end) = marktab(1, fillmarkpos+1:end) + filltime;
        end
        mask = all(abs(diff(pulsetab(2:end, :), [], 2)) < 1e-14);
        pulsetab(:, [false, mask(2:end) & mask(1:end-1)]) = [];
        pulse.data = struct;
        pulse.data.pulsetab = pulsetab;
        pulse.data.marktab = sortrows(marktab',1)';
        pulse.data.readout = readout;
        pulse.data.elem=pulseDef;  % Copy forward documentation.
        pulse.format = 'tab';        
    otherwise
        error('Invalid format %s.', pulse.format);
end
end  