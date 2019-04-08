function [pars, chisq, cov,exitflag,fitinfo] = mfitwrap(data,model,beta0,opts,mask)
% This function fits multiple datasets simultaneously using non-linear least squares
% function [ p, chisq, cov ] = mfitwrap(data,model,beta0,opts,mask)
% data is a struct array of length ndatasets with fields
%  data.x  -- x coordinates
%  data.y  -- y coordinates
%  data.vary -- y variances.  Assumed 1 if not provided.
%  data.rng -- range to fit in x (optional, empty means everything)
% model is a struct array with length ndatasets of models to use in the fitting.
%  model.fn  -- fit function, of the form fn(p,x) where p are the parameters
%  model.yfn -- function that transforms y->yfn(p,y). (e.g.if you want to take log of data)
%  model.pt  -- parameter trafofn. (so if 5th model uses params 1 and 7 in
%  mfitwrap as parameters 1 and 2 in its fit function, model(5).pt = @(p) [p(1) p(7)];
% mask - set to 1 to allow parameter to vary, 0 to hold fixed. Default is all 1.
% pars is a list of fitted parameters
% chisq is the reduced chisq; chi^2/(npts-ndof)
% cov is the covariance of the fit parameters
% exitflag : reason solver stopped.
% fitInfo : information about fit process.
% opts can include:
%   plinit : plot initial guess
%   plfit : plot fit
%   optimplot:
%   lm: use levenberg-marquardt
%   fine: use lower tolerances
%   samefig: plot on gcf
%   nofit: don't fit data.
% default: plfit plinit optimplot
% add control of MaxIterations, MaxFunctionEvaluations
if ~exist('opts','var') || isempty(opts), opts='plfit plinit optimplot'; end
if ~exist('mask','var'), mask=true(size(beta0)); end
mask = logical(mask);
fitOpts=optimoptions(optimoptions('lsqnonlin'),'display','off'); % could use 'iter' for display to 'display output at each iteration, and give the default exit message'
%FIxme:  Is this default? 
if isopt(opts,'optimplot'), fitOpts=optimoptions(fitOpts,'PlotFcns',@optimplotresnorm); end % plots the norm of the residuals
if isopt(opts,'lm'), fitOpts=optimoptions(fitOpts,'Algorithm','levenberg-marquardt'); end
fitOpts = optimoptions(fitOpts,'MaxIterations',300,'MaxFunctionEvaluations',10000);

if isopt(opts,'fine'), fitOpts=optimoptions(fitOpts,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e5,'MaxIter',1e4); end

for i =1:length(data) % Grab y variance.
    if ~isfield(data(i),'vary') || isempty(data(i).vary) 
        data(i).vary=ones(size(data(i).y));
    end
    data(i).stdy = sqrt(data(i).vary); 
end

if isopt(opts,'plinit') % Plot initial guess
    f=figure(60); clf;
    f.Name='Initial Guess';
    lsqfunFull(data,model,beta0,['mustplot samefig ' opts],beta0,true(size(beta0)));        
end 
if ~isopt(opts,'nofit') % Perform fit
    % lsqnonlin tries to minimize function given.
    if (isfield(model,'pt') && ~isempty(model(i).pt)) || (isfield(model,'yfn') && ~isempty(model(i).yfn))
        [fitPars, chisq,~, exitflag, fitinfo,~, jac] = lsqnonlin(@(p) lsqfunFull(data,model,p,opts,beta0,mask), beta0(mask),[],[],fitOpts);
    else
        [fitPars, chisq,~, exitflag, fitinfo,~, jac] = lsqnonlin(@(p) lsqfun(data,model,p,beta0,mask), beta0(mask),[],[],fitOpts);
    end
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
if isopt(opts,'plfit') && ~isopt(opts,'nofit') % Plot best fit.
    f=figure(61); clf;
    f.Name = 'Best Fit';
    lsqfunFull(data,model,pars,['mustplot samefig ' opts],pars,true(size(pars)));    
end
end


function err=lsqfun(data, model, fitPars,beta0,mask)
% Returns residual for current set of fitPars.
pars(mask)=fitPars; % Pars set to current best fit value.
pars(~mask)=beta0(~mask); % Unfit paramaters are set to initial guess.
err = [];

for i=1:length(data)    
    fitData=model(i).fn(pars, data(i).x); % Get current value of model.
    err = [err (fitData-data(i).y)./data(i).stdy]; %#ok<*AGROW> Add error divided by variance to list.    
end
err=err';
end

function err=lsqfunFull(data, model, fitPars,opts,beta0,mask)
% Returns residual for current set of fitPars.
persistent lastplot; %persistent num
pars(mask)=fitPars; % Pars set to current best fit value.
pars(~mask)=beta0(~mask); % Unfit paramaters are set to initial guess.
err = [];

doplot = isopt(opts,'mustplot');
if isopt(opts,'plotiter') &&( isempty(lastplot) || (now > lastplot + 15/(24*60*60)))
    doplot = 1;
end
if doplot
    lastplot = now;
    if ~isopt(opts,'samefig')
        f=figure(62);
        if ~isopt(opts,'noclear'), clf; end
        f.Name = 'Iteration Display';
    end
    rows=max(floor(sqrt(length(data))*.85),1);
    cols=ceil(length(data)/rows);
    %fprintf('%d \n',num);    
    ga = tightSubplot([rows,cols]);    
end

for i=1:length(data)
    if isfield(model,'pt') && ~isempty(model(i).pt) % Grab parameters for current model.
        parsCurrModel=model(i).pt(pars);
    else
        parsCurrModel=pars;
    end
    
    if isfield(model(i),'yfn') % Get y data.
        [y,vary]=model(i).yfn(pars,data(i).y,vary);
        if any(vary < 0)
            error('Negative variance');
        end
    else
        y=data(i).y;
    end
    fitData=model(i).fn(parsCurrModel, data(i).x); % Get current value of model.
    err = [err (fitData-y)./data(i).stdy]; %#ok<*AGROW> Add error divided by variance to list.
    if doplot        
        if isopt(opts,'err')
            if isopt(opts,'green')
                errorbar(ga(i),data(i).x,y,data(i).stdy,'g');
            else
                errorbar(ga(i),data(i).x,y,data(i).stdy,'.-','CapSize',1);
            end
        else
            plot(ga(i),data(i).x,y,'k.-');
        end
        hold(ga(i),'on');
        
        if any(imag(fitData) ~= 0), error('Imaginary fit');   end
        if ~isopt(opts,'nofunc'), plot(ga(i),data(i).x,fitData);  end
    end
end
err=err';
end