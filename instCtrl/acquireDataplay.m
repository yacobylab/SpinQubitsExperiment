function [result] = acquireDataplay(boardHandle,npoints,rate)
% Make an AutoDMA acquisition from dual-ported memory.
% global variable set in configureBoard.m
global SamplesPerSec


result = false; % set default return code to indicate failure
AlazarDefs %call mfile with library definitions
acquisitionLength_sec = npoints/rate; % Select the total acquisition length in seconds (or 0 to acquire until aborted)
samplesPerChannel = 1024000; % Select the number of samples per channel in each buffer
bufferTimeout_ms = 5000; % This is the amount of time to wait for for each buffer to be filled
saveData = true; % Select if you wish to save the sample data to a binary file
drawData = false; % TODO: Select if you wish to plot the data to a chart
channelMask = CHANNEL_A + CHANNEL_B; % TODO: Select which channels read from on-board memory (A, B, C, D, or all)

% Calculate the number of enabled channels from the channel mask 
channelCount = 0;
channelsPerBoard = 2;
for channel = 0:channelsPerBoard - 1
    channelId = 2^channel;
    if bitand(channelId, channelMask)
        channelCount = channelCount + 1;
    end
end
if (channelCount < 1) || (channelCount > channelsPerBoard)
    fprintf('Error: Invalid channel mask %08X\n', channelMask);
    return
end


[retCode, boardHandle, maxSamplesPerRecord, bitsPerSample] = daqfn('AlazarGetChannelInfo', boardHandle, 0, 0); % Get the sample and memory size
if retCode ~= ApiSuccess
    fprintf('Error: AlazarGetChannelInfo failed -- %s\n', errorToText(retCode));
    return
end
if samplesPerChannel > maxSamplesPerRecord
    fprintf('Error: Too many samples per buffer %u max %u\n', samplesPerChannel, maxSamplesPerRecord);
    return
end

% Calculate the size of each buffer in bytes
bytesPerSample = floor((double(bitsPerSample) + 7) / double(8));
samplesPerBuffer = samplesPerChannel * channelCount;
bytesPerBuffer = uint32(bytesPerSample) * samplesPerBuffer;

% Find the number of buffers in the acquisition
if acquisitionLength_sec > 0 
    samplesPerAcquisition = uint32(floor((SamplesPerSec * acquisitionLength_sec + 0.5)));
    buffersPerAcquisition = uint32(floor((samplesPerAcquisition + samplesPerChannel - 1) / samplesPerChannel));
else
    buffersPerAcquisition = hex2dec('7FFFFFFF');  % acquire until aborted
end

% TODO: Select the number of DMA buffers to allocate.
% The number of DMA buffers must be greater than 2 to allow a board to DMA into
% one buffer while, at the same time, your application processes another buffer.
bufferCount = uint32(8);

% Create an array of DMA buffers 
buffers = cell(1, bufferCount);
for j = 1 : bufferCount
    pbuffer = calllib('ATSApi', 'AlazarAllocBufferU16', boardHandle, samplesPerBuffer);
    if pbuffer == 0
        fprintf('Error: AlazarAllocBufferU16 %u samples failed\n', samplesPerBuffer);
        return
    end    
    buffers(1, j) = { pbuffer };
end

% TODO: Select AutoDMA flags as required
% ADMA_EXTERNAL_STARTCAPTURE - call AlazarStartCapture to begin the acquisition
% ADMA_TRIGGERED_STREAMING - acquire a single gapless record spanning multiple buffers
%   after a trigger event.
% ADMA_INTERLEAVE_SAMPLES - interleave samples for highest throughput
admaFlags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRIGGERED_STREAMING + ADMA_INTERLEAVE_SAMPLES;
% Configure the board to make an AutoDMA acquisition
retCode = calllib('ATSApi', 'AlazarBeforeAsyncRead', boardHandle, channelMask, 0, samplesPerChannel, 1, buffersPerAcquisition, admaFlags);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(retCode));
    return
end

% Post the buffers to the board
for bufferIndex = 1 : bufferCount
    pbuffer = buffers{1, bufferIndex};
    retCode = calllib('ATSApi', 'AlazarPostAsyncBuffer', boardHandle, pbuffer, bytesPerBuffer);
    if retCode ~= ApiSuccess
        fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode));
        return
    end        
end

% Update status
if buffersPerAcquisition == hex2dec('7FFFFFFF')
    fprintf('Capturing buffers until aborted...\n');
else
    fprintf('Capturing %u buffers ...\n', buffersPerAcquisition);
end

% Arm the board system to begin the acquisition 
retCode = calllib('ATSApi', 'AlazarStartCapture', boardHandle);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarStartCapture failed -- %s\n', errorToText(retCode));
    return;
end

% Wait for sufficient data to arrive to fill a buffer, process the buffer,
% and repeat until the acquisition is complete
startTickCount = tic;
updateTickCount = tic;
updateInterval_sec = 0.1;
buffersCompleted = 0;
captureDone = false;
success = false;

while ~captureDone
    
    bufferIndex = mod(buffersCompleted, bufferCount) + 1;
    pbuffer = buffers{1, bufferIndex};
    
    % Wait for the first available buffer to be filled by the board
    [retCode, boardHandle, bufferOut] = ...
        calllib('ATSApi', 'AlazarWaitAsyncBufferComplete', boardHandle, pbuffer, bufferTimeout_ms);
    if retCode == ApiSuccess
        % This buffer is full
        bufferFull = true;
        captureDone = false;
    elseif retCode == ApiWaitTimeout
        % The wait timeout expired before this buffer was filled.
        % The board may not be triggering, or the timeout period may be too short.
        fprintf('Error: AlazarWaitAsyncBufferComplete timeout -- Verify trigger!\n');
        bufferFull = false;
        captureDone = true;
    else
        % The acquisition failed
        fprintf('Error: AlazarWaitAsyncBufferComplete failed -- %s\n', errorToText(retCode));
        bufferFull = false;
        captureDone = true;
    end
    
    if bufferFull
        % TODO: Process sample data in this buffer.
        %
        % NOTE:
        % While you are processing this buffer, the board is already
        % filling the next available DMA buffer.
        %
        % You must finish processing this buffer before the board fills
        % all of its available DMA buffers and on-board memory.
        %
        % Records are arranged in the buffer as follows: R0A, R0B
        %
        % Samples values are arranged contiguously in each record.
        % A 14-bit sample code is stored in the most significant bits of
        % each 16-bit sample value.
        %
        % Sample codes are unsigned by default where:
        % - 0x0000 represents a negative full scale input signal;
        % - 0x2000 represents a 0V signal;
        % - 0x3fff represents a positive full scale input signal.
        
        setdatatype(bufferOut, 'uint16Ptr', 1, samplesPerBuffer);
        
        % Save the buffer to file
        if fid ~= -1
            samplesWritten = fwrite(fid, bufferOut.Value, 'uint16');
            if samplesWritten ~= samplesPerBuffer
                fprintf('Error: Write buffer %u failed\n', buffersCompleted);
            end
        end
        
        % Display the buffer on screen
        if drawData
            plot(bufferOut.Value);
        end
        
        % Make the buffer available to be filled again by the board
        retCode = calllib('ATSApi', 'AlazarPostAsyncBuffer', boardHandle, pbuffer, bytesPerBuffer);
        if retCode ~= ApiSuccess
            fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(retCode));
            captureDone = true;
        end
        
        % Update progress
        buffersCompleted = buffersCompleted + 1;
        if buffersCompleted >= buffersPerAcquisition
            captureDone = true;
            success = true;
        elseif toc(updateTickCount) > updateInterval_sec
            updateTickCount = tic;
            
        end
        
    end % if bufferFull    
end

    fprintf('Captured %u buffers in %g sec (%g buffers per sec)\n', buffersCompleted, transferTime_sec, buffersPerSec);
    fprintf('Transferred %u bytes (%.4g  per sec)\n', bytesTransferred, bytesPerSec);   

% Save the transfer time
transferTime_sec = toc(startTickCount);

% Close progress window
delete(waitbarHandle);

% Abort the acquisition
retCode = calllib('ATSApi', 'AlazarAbortAsyncRead', boardHandle);
if retCode ~= ApiSuccess
    fprintf('Error: AlazarAbortAsyncRead failed -- %s\n', errorToText(retCode));
end

% Close the data file
if fid ~= -1
    fclose(fid);
end

% Release the buffers
for bufferIndex = 1:bufferCount
    pbuffer = buffers{1, bufferIndex};
    retCode = calllib('ATSApi', 'AlazarFreeBufferU16', boardHandle, pbuffer);
    if retCode ~= ApiSuccess
        fprintf('Error: AlazarFreeBufferU16 failed -- %s\n', errorToText(retCode));
    end
    clear pbuffer;
end

% Display results
if buffersCompleted > 0 
    bytesTransferred = double(buffersCompleted) * double(bytesPerBuffer);
    if transferTime_sec > 0 
        buffersPerSec = buffersCompleted / transferTime_sec;
        bytesPerSec = bytesTransferred / transferTime_sec;
    else
        buffersPerSec = 0;
        bytesPerSec = 0;
    end
end
% set return code to indicate success
result = success;

end