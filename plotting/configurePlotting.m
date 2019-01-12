%% script that controls the autoplotting of a folder. 
global pptdata; 

%otherCat{1} = {'NT', 'Wall', 'Wall4k', 'chrg', 'chrg_sens', 'scanZoom', 'sensLock', 'sensorD', 'sensorDAQ', 'DAQsensor', };
%otherCat{2} = {'sit', '3a', '3b' ,'3b_T34', '4a', '4b', '4b_T34', 'N34', 'SD4bot', 'SD4mid', 'SD4top', 'line_2a','line_SD4top'};

dev(1).IndRng = [1 339; 391 478; 479 799; 800 845; 846 1438]; dev(1).name ='res';
%dev(1).CooldownDesc = {'F bad Amp'; 'F bad trap'; 'F ground'; 'F short'; 'F tune'}; 
dev(1).shtName = 'res'; 
dev(1).CooldownDesc = {'F Amp'; 'F Trap'; 'F Grnd'; 'F Blow'; 'F Tune'}; 

dev(2).IndRng = [1439 1596; 2745 2801]; dev(2).name = 'Ozone';%Oz 
dev(2).CooldownDesc = {'F Tune'; '4 Hyst'};
dev(2).shtName = 'oz'; 

dev(3).IndRng = [1597 1612; 1613 1694; 1695 1803; 2116 2247; 2300 2408; 2709 2744]; dev(3).name = 'EBeam';% Ebeam 
dev(3).CooldownDesc = {'F missGt'; '4 testD1';'4 testD2'; 'F tune'; 'F tune2'; '4 hyst'};
dev(3).shtName = 'EB'; 

dev(4).IndRng = [1806 1874; 1891 1936; 1937 2115]; dev(4).name = 'Old'; % Old 
dev(4).CooldownDesc = {'4 test'; 'F noisy'; 'F tune'}; 
dev(4).shtName = 'old'; 

dev(5).IndRng = [2409 2428; 2429 2455; 2456 2708; 2802 2837]; dev(5).qpcRng = [86 145]; dev(5).name = 'EBeamRes1';% EBeam Res 1 
dev(5).CooldownDesc = {'4 test1'; '4 test2'; 'F tune'; 'F hyst'}; 
dev(5).qpcDesc = {'F hyst'}; dev(5).shtName = 'EBR1'; 

dev(6).IndRng = [2844 2859]; dev(6).qpcRng = [1 34; 35 44; 45 85]; dev(6).name = '3qubit';  % 3 qubit 
dev(6).CooldownDesc = {'4 hyst'}; dev(6).qpcDesc = {'4 hyst','4 high','4 hyst'};
dev(6).shtName = '3Q'; 


dev(7).IndRng = [2861 2874; 2877 2882; 2887 2894; 2895 2920; 2924 NaN]; 
dev(7).qpcRng = [147 158; 162 164; 165 183; 184 191; 192 226; 227 241; 242 275; 276 NaN]; dev(7).name = 'EBeamRes2'; % EBeam Res 2 
dev(7).CooldownDesc = {'4 test', '4 test2','4 test100','4test1002','F tune'};
dev(7).qpcDesc = {'4 a', '4 ohmc', '4 unc', '4 chkGt','F4 hystB','F42 hystB','FC hystA', 'F4 hystA'}; 
dev(7).shtName = 'EBR2'; 


dev(8).IndRng = [2248 2298]; dev(8).name = 'Noisy';%noisy, weird. 
dev(8).CooldownDesc = {'4'}; dev(8).shtName = 'Nois'; 
pptdata.dev = dev; 

cats{1} = {'qpc_','midSweep%s','line_.*','sit%s',}; cats{2} = {'hyst%s.*','qpc%s','qpc_4k_%s','oneGate%s','singGate%s'};
cats{3} = {'NT%s','QPCone%s'}; 
cats{4} = {'twoD%s.*','offT%s.*','Wall%s','Wall4k%s','scanWall%s','Leads%s','scanZoom','singDot%s','SDLock%s','SD%sLock','scanSD%s','SD%s','sensDot%s','sensorLock%s','sensor%s','sensorD%s','sensorDAQ%s','DAQsensor%s'};
cats{5} = {'sensLock%s','sens%s','sens_%s','sens%sRot','chrg_sens','junc%s','chrg%s'}; 
cats{6} = {'RF%s'}; 
scanList = cats{:}; 
pptdata.cats = cats;
pptdata.dataFolder = 'Z:\qDots\data\data_2015_11_05'; 
pptdata.qpcFolder = 'qpc_2016_05_16\'; 
pptdata.folder = 'Z:\qDots\PPT\ppt_2015_11_05\';

n=9;
dev(n).IndRng = [3295 3307; 3308 3390; 3391 3474; 3475, 3929; 3930 3932; 3933 4317]; 
dev(n).qpcRng = [318 338; 339 367; 369 376; 377 447; 448 462; 463 524; 525 532; 533 557]; 
dev(n).name = 'SH1-d-M5-22-13.1'; 
dev(n).CooldownDesc = {'init4','fridgeNoBias','fridge','fridgeSensBias','fridge200','fridge235'}; 
dev(n).qpcDesc = {'init4','RT','Umansky','fridgeNoBias','fridge100','fridgeSensBias100','fridge200','fridge235'}; 
dev(n).shtName = 'F1251';
pptdata.dev(9)=dev(9); 
pptdata.filename = sprintf('dev%d',n); 


n = 10;
dev(n).IndRng = [4318 6318];
dev(n).qpcRng = [558 602; 603 644];
dev(n).name = 'SH1-d-M5-22-13.2';
dev(n).CooldownDesc = {'initFridge'};
dev(n).qpcDesc = {'4k','fridge'};
dev(n).shtName = 'F1252';
pptdata.dev(n)=dev(n); 

n = 11; 
dev(n).IndRng = [6319 7156];
dev(n).qpcRng = [645 698; 699 722];
dev(n).name = 'CB1-c-M5-22-13.1';
dev(n).CooldownDesc = {'initFridge'};
dev(n).qpcDesc = {'4k','fridge'};
dev(n).shtName = 'ALD1';
pptdata.dev(n)=dev(n); 

n = 12; 
dev(n).IndRng = [7157 7168];
dev(n).qpcRng = [723 873];
dev(n).name = 'CB2-c-M5-22-13.1';
dev(n).CooldownDesc = {'initFridge'};
dev(n).qpcDesc = {'4k','fridge'};
dev(n).shtName = 'ALD2';
pptdata.dev(n)=dev(n); 

n = 13; 
dev(n).IndRng = [7169 7294; 7295 NaN];
dev(n).qpcRng = [874 930; 931 940; 941 NaN];
dev(n).name = 'CB2-d-M5-22-13.1';
dev(n).CooldownDesc = {'initFridge','fridge2'};
dev(n).qpcDesc = {'4k','fridge','fridge2'};
dev(n).shtName = 'ALD3';
pptdata.dev(n)=dev(n); 
save('Z:/Shannon/Data/pptdata','pptdata'); 
%% Stuff to update each time matlab starts 

pptdata.dir = dir; 
pptdata.qpcDir = dir(pptdata.qpcFolder); 

% should we be using that one more often???
pptdata.fileNames = sortFiles; % Cell array of names in order by datenum.
pptdata.qpcfileNames = sortFiles(pptdata.qpcFolder); % Cell array of names in order by datenum.
save('Z:/Shannon/Data/pptdata','pptdata'); 
%% we want to plot all of the qpc stuff everything else.
out={};
scanTypes = length(scanList);
%pptControl('start')
for i = 1:length(dev)
    for j = 1:size(dev(i).qpcRng,1)
        out{end+1}=autoPlotFunc(i,-j,1,'hyst old');
    end
    for j = 1:size(dev(i).IndRng,1)
        out{end+1}=autoPlotFunc(i,j,1,'hyst old');
    end    
end
%%
out={};
pptControl('start')
for i = 1:length(dev)
    for j = 1:size(dev(i).qpcRng,1) 
        out{end+1}=autoPlotFunc(i,-j,5,'old');
    end
    for j = 1:size(dev(i).IndRng,1)
        out{end+1}=autoPlotFunc(i,j,5,'old');
    end    
end
%pptdata.filename = 'qpcHyst2';
%pptControl('save'); 
%%
for i =1:1
    pptControl('start')
    for j = 1:3
        for k = 1:size(dev(i).IndRng,1)
            autoPlotFunc(i,k,'',j,'old color ');
            pptdata.filename = sprintf('dev%d0',i);            
        end
    end
    pptControl('save'); pptControl('end');
end
