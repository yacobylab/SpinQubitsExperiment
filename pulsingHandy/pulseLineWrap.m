function pulseLineWrap(scan,pulse)
scan.consts(end+1).setchan = 'PulseLine'; 
scan.consts(end).val = awgseqind(pulse)+1; 
end