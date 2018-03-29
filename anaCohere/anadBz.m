function anadBz
figure(17); clf; hold on;
a = procPlsData; 

beta0 = [.5 .4 .1 0.03*2*pi 200];
fitfn = @(y,x) y(1)+y(2)*sin(y(4)*x+y(3)).*exp(-x.^2 ./ y(5).^2);
data = squeeze(nanmean(a.data{1}));
% offset, amplitude, freq (rad/s), phase, t2*)
pars = fitwrap('plfit',a.xv{1},data',beta0,fitfn);
fprintf('T2*: %3.1f \n',pars(5)); 
end