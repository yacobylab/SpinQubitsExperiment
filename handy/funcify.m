function func=funcify(fitFn)
if ischar(fitFn)
    func = str2func(fitFn);
else
    func = fitFn;
end
end