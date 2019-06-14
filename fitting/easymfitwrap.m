function easymfitwrap(opts,x,y,beta0,fitfn,mask)
data.x = x; data.y = y; 
model.fn = fitfn; 

if isstruct(beta0) && ~isempty(beta0.fn)
    beta0 = beta0.fn(x, y, beta0.args{:});
end
if ischar(model.fn), model.fn=str2func(model.fn); end

mfitwrap(data,model,beta0,opts,mask);
end