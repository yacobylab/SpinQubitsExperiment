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
%   err - plot error bars
%   green 
%   plinit : plot initial guess
%   plfit : plot fit 
%   optimplot
%   lm: use levenberg-marquardt
%   fine: use lower tolerances
%   samefig: plot on gcf
%   mustplot
%   nofit
%   plotiter
%   noclear
%   robust
%   nofunc
% default: plfit plinit optimplot
global fbdata 
%if ~exist('opts','var') || isempty(opts), opts='plfit plinit optimplot'; end
if ~exist('mask','var'), mask=true(size(beta0)); end
mask = logical(mask); 
fitOpts=fbdata.fitOpts; 
if isopt(opts,'optimplot'), fitOpts=optimset(fitOpts,'PlotFcns',@optimplotresnorm); end
if isopt(opts,'lm'), fitOpts=optimset(fitOpts,'Algorithm','levenberg-marquardt'); end
if isopt(opts,'fine'), fitOpts=optimset(fitOpts,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e5,'MaxIter',1e4); end
if isopt(opts,'plinit')
    f=figure(60); clf;
    f.Name='Initial Guess';
    lsqfun(data,model,beta0,['mustplot samefig' opts],beta0,true(size(beta0)));
    
    f=figure(63); clf;
    f.Name='Comparison';
    lsqfun(data,model,beta0,['mustplot samefig nofunc' opts],beta0,true(size(beta0)));
end
if ~isopt(opts,'nofit')
    [fitPars, chisq, ~, ~,~, ~, jac] = lsqnonlin(@(p) lsqfun(data,model,p,opts,beta0,mask), beta0(mask),[],[],fitOpts);
    covt = pinv(full(jac' * jac));  % Should this be inv not pinv?  singularity implies some fit paramteres don't matter....
    cov=zeros(length(beta0),length(beta0));
    cov(find(mask),find(mask))=covt; %#ok<FNDSB>
    pars=beta0;
    pars(mask)=fitPars;
    npts=numel([data.x]);
    chisq=chisq/(npts-sum(mask));
else
    pars=beta0;
    chisq=0;
end
if isopt(opts,'plfit') && ~isopt(opts,'nofit')
    f=figure(61); clf;
    f.Name = 'Best Fit';
    lsqfun(data,model,pars,['mustplot samefig' opts],pars,true(size(pars)));
    
    figure(63);
    lsqfun(data,model,beta0,['mustplot samefig noclear nofunc' opts],beta0,true(size(beta0)));
end
end

function err=lsqfun(data, model, fitPars,opts,beta0,mask)
persistent lastplot;
pars(mask)=fitPars;
pars(~mask)=beta0(~mask);
err = [];
doplot = isopt(opts,'mustplot');
if isopt(opts,'plotiter') &&( isempty(lastplot) || (now > lastplot + 0.5/(24*60*60)))
    doplot = 1;
end
if doplot
    lastplot = now;
    if ~isopt(opts,'samefig')
        f=figure(62);
        if ~isopt(opts,'noclear')
            clf;
        end
        f.Name = 'Iteration Display';
    end
    rows=floor(sqrt(length(data)));
    cols=ceil(length(data)/rows);
end
for i=1:length(data)
    if isfield(model,'pt') && ~isempty(model(i).pt)
        parsCurrModel=model(i).pt(pars);
    else
        parsCurrModel=pars;
    end
    if ~isfield(data(i),'vary') || isempty(data(i).vary)
        sy=ones(size(data(i).y));
    else
        sy=data(i).vary;
    end
    if isfield(model(i),'yfn')
        [y,sy]=model(i).yfn(pars,data(i).y,sy);
        if any(sy < 0)
            error('Negative variance');
        end
    else
        y=data(i).y;
    end
    fitData=model(i).fn(parsCurrModel, data(i).x);
    err = [err (fitData-y)./sqrt(sy)];
    if any(imag(err) ~= 0)
        error('Imaginary error');
    end
    if doplot
        subplot(rows,cols,i);
        if isopt(opts,'err')
            if isopt(opts,'green')
                errorbar(data(i).x,y,sqrt(sy),'g');
            else
                errorbar(data(i).x,y,sqrt(sy),'r');
            end
        end
        hold on;
        plot(data(i).x,y,'kx-');                       
        if any(imag(fitData) ~= 0), error('Imaginary fit');   end
        if ~isopt(opts,'nofunc'), plot(data(i).x,fitData,'b-');  end
    end
end
if doplot, drawnow; end
err=err';
if isopt(opts,'robust'), err = err ./ sqrt(abs(err)); end % robust fit. 
end