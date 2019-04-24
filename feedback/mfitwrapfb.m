function [pars, chisq, cov] = mfitwrapfb(data,model,beta0,opts,mask)
% limited fast version of mfitwrap to use for fb. 
% function [ p, chisq, cov ] = mfitwrapfb(data,model,beta0,opts,mask)
% This function fits multiple datasets simultaneously using lsqnonlin
% Plots : 
% data is a struct array of length ndatasets with fields
%  data.x  -- x coordinates
%  data.y  -- y coordinates
%  data.vary -- y variances.  Assumed 1 if not provided.
%  data.rng -- range to fit in x (optional, empty means everything)
% model is a struct array with length ndatasets of models to use in the fitting.
%  model.fn  -- fit function, of the form fn(p,x) where p are the parameters
%  model.yfn -- bizarre fit function that transforms y->yfn(p,y). (e.g.if you want to take log of data)
%  model.pt  -- parameter trafofn. (so if 5th model uses params 1 and 7, model(5) = @(p) [p(1) p(7)]; 
% mask - set to 1 to allow parameter to vary, 0 to hold fixed.
% p is the output parameters
% chisq is the reduced chisq; chi^2/(npts-ndof)
% cov is the covariance of the fit parameters
% opts can include:
%   lm: use levenberg-marquardt
%   fine: use lower tolerances
%   robust
global fbdata 
if ~exist('mask','var'), mask=true(size(beta0)); end
mask = logical(mask); 
fitOpts=fbdata.fitOpts; 

if isopt(opts,'lm'), fitOpts=optimset(fitOpts,'Algorithm','levenberg-marquardt'); end
if isopt(opts,'fine'), fitOpts=optimset(fitOpts,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e5,'MaxIter',1e4); end

for i =1:length(data) % Grab y variance.
    if ~isfield(data(i),'vary') || isempty(data(i).vary) 
        data(i).vary=ones(size(data(i).y));
    end
    data(i).stdy = sqrt(data(i).vary); 
end

[fitPars, chisq, ~, ~,~, ~, jac] = lsqnonlin(@(p) lsqfun(data,model,p,beta0,mask), beta0(mask),[],[],fitOpts);
covt = pinv(full(jac' * jac));  % Should this be inv not pinv?  singularity implies some fit paramteres don't matter....
cov=zeros(length(beta0),length(beta0));
cov(find(mask),find(mask))=covt; %#ok<FNDSB>
pars=beta0;
pars(mask)=fitPars;
npts=numel([data.x]);
chisq=chisq/(npts-sum(mask));
end

function err=lsqfun(data, model, fitPars,beta0,mask)
pars(mask)=fitPars;
pars(~mask)=beta0(~mask);
err = [];

for i=1:length(data)
    fitData=model(i).fn(pars, data(i).x);
    err = [err (fitData-data(i).y)./data(i).stdy];    
end
err=err';
end