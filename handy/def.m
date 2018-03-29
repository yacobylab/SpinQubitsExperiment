function strct=def(strct,fldname,val)
if ~isfield(strct,fldname)
    strct.(fldname)=val;
end
end