function clearBuffers
global smdata; 
ico = inl('ATS660');

missedbuf = [];
for j = 1:length(smdata.inst(ico(1)).data.buffers) % Free buffers
    try
        daqfn('FreeBufferU16', boardHandle, smdata.inst(ico(1)).data.buffers{j});
    catch
        missedbuf(end+1)=j; %#ok<AGROW>
    end
end

end