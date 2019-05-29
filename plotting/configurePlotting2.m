n = 1; 
dev(n).IndRng = [1 NaN];
dev(n).qpcRng = [1 NaN];
dev(n).tuneRng = [1 NaN];
dev(n).name = 'CB1-M07-02-13.1';
dev(n).CooldownDesc = {'fridge1'};
dev(n).qpcDesc = {'4k','fridge1'};
dev(n).shtName = 'HSQ1';
pptdata.dev = dev; 

cats{1} = {'qpc_','midSweep%s','line_.*','sit%s',}; cats{2} = {'hyst%s.*','qpc%s','qpc_4k_%s','oneGate%s','singGate%s'};
cats{3} = {'NT%s','QPCone%s'}; 
cats{4} = {'twoD%s.*','offT%s.*','Wall%s','Wall4k%s','scanWall%s','Leads%s','scanZoom','singDot%s','SDLock%s','SD%sLock','scanSD%s','SD%s','sensDot%s','sensorLock%s','sensor%s','sensorD%s','sensorDAQ%s','DAQsensor%s'};
cats{5} = {'sensLock%s','sens%s','sens_%s','sens%sRot','chrg_sens','junc%s','chrg%s'}; 
cats{6} = {'RF%s'}; 
scanList = cats{:}; 
pptdata.cats = cats;

%%
createDevppt('R')