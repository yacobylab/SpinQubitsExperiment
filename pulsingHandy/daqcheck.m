% get more alazar info
inst = sminstlookup('ATS660'); 
boardh = smdata.inst(inst).data.handle; 
param_vals=[268435513,268435520, 268435522,268435524,268435525,268435536,268435537,268435538]; 
pparam=libpointer('int32Ptr',int32(1)); 
ch=uint8(1); 
msg={'Size in bytes','Number of allocated buffers', 'Data signed','Samples per timestep'...
  'Number records captured','Number buffers queued','Number buffers full','Number buffers empty'}; 
for i = 1:length(param_vals)
param=uint32(param_vals(i)); 
 %  daqfn('GetParameter',smdata.inst(11).data.handle,ch,param,pparam)
 try
   daqfn('GetParameter',boardh,ch,param,pparam)
   fprintf('%s : %d \n', msg{i}, pparam.value)       
 end
end
        % 1 / data width: 14 
        % 1/ bites in each dma : 0 
        % # of allocated buffers: 64 [none in data] 
        % # of sample clocks per timestep clock: 2 
        % 268435561: size in bytes of each DMA 
        % 268435566 : get # allocated buffers
        % 268435568 : data format (signed/unsigned) 
        % #s: 268435570, samples per timestep 
        % 268435571 : # records captured since start of acq
        % 268435581 : # buffers queued 
        % 268435582 : # buffers full 
        % 268435583 : # buffers empty 
        % 
        
     maxrec=uint32(1); 
     pmaxrec=libpointer('uint32Ptr',maxrec); 
     daqfn('GetMaxRecordsCapable', boardh,1024,pmaxrec);
   fprintf('Max Records: %d \n',maxrec)
   
 pmaj=libpointer('uint8Ptr',uint8(1));
 pmin=libpointer('uint8Ptr',uint8(1));
 prev=libpointer('uint8Ptr',uint8(1));
 daqfn('GetSDKVersion',pmaj,pmin,prev);
 
 fprintf('Driver version is %d.%d.%d \n', pmaj.value, pmin.value, prev.value);
 nboard = daqfnAns('BoardsFound'); 
 
 pmaj=libpointer('uint32Ptr',uint8(1));
 pmin=libpointer('uint8Ptr',uint8(1));
 daqfn('GetChannelInfo',boardh,pmaj,pmin);
 fprintf('Memory Size %d, nbits %d \n',pmaj.value,pmin.value); 
  a = daqfnAns('GetStatus',boardh)
%   0x00000001 - At least 1 trigger occurred.
% 0x00000002 - Channel A input went beyond the limits of its input
% range for at least one sample during the acquisition.
% 0x00000004 - Channel B input went beyond the limits of its input
% range for at least one sample during the acquisition
% 0x00000008 – PLL locked (ATS660 and ATS9462 Only)
% 0xFFFFFFFF – Error – No value is valid.
        %%
t = daqfnAns('Triggered',boardh)
        
t = daqfnAns('GetWhoTriggeredBySystemHandle',boardh,1,1)
        %Gather Alazar ifno. 
try
    val=daqfn('GetBoardKind', boardh); %didn't work Err 16
end
     %no, it said it was a 9440. 
 %% didn't work -- failed in callib. 
 
         % bits per sample is 14 
        %memory size is 1.34217728e8 
 %%   
    
      
        % Get the SDK version 
        % 268435570 Samples clocks per time step 
        
        % 5.9.0 - version of driver library file, not SDK version number. 
        
 %% System Overview 
 
 %nboard=daqfn_ans('BoardsFound')
 %typboard=daqfn_ans('GetBoardKind',smdata.inst(8).data.handle);
 %fprintf('%d Boards found of type %d \n',nboard,typboard)
%   
%  pmemsize=libpointer('uint32Ptr',uint32(zeros(1,1)));
%  pbitssamp=libpointer('uint8Ptr',uint8(zeros(1,1)));
%  daqfn('SetParameter',smdata.inst(11).data.handle,0,hex2dec('10000040'),32)
