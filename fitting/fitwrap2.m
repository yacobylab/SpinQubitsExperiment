function [beta1,res,jac,COVB,mse,err] = fitwrap(opts, x, y, beta0, model, mask)
% Performs nlinfit with plotting, a few additional options. 
% function [beta1,r,j,COVB,mse,err] = fitwrap(ctrl, x, y, beta0, model, mask)
% model defined as  y = model(beta, x), can be a string or function handle.
% ctrl: plinit, plfit, woff, nofit, pause, samefig, fine, robust
% ctrls:
%   plinit: plots the initial guess
%   plfit: plots the fit
%   woff: turns off warnings
%   nofit: does not fit the data
%   pause: will pause after each plot (until press any key)
%   samefig: uses current figure. Otherwise clears 500 and uses it. 
%   fine: fitting option
%   robust: fitting option
%   resid: Plot residual. 
% beta0 can be a vector, cell array with initialization function handles
% or initial value vectors, or a struct array with fields fn and args, and optionally vals.
% Finite vals entries override the return of the initialization function.
% Initialization functions are called as beta0 = fn(x, y, args).
% data should be in columns.
% Outputs:
%   beta1 are the fitted parameters
%   res are the residuals
%   jac is the Jacobian
%   COVB is the variance-covariance matrix
%   mse is the mean square error. 
%   err are the upper and lower error bounds for the predictions. the format is val=beta(i,j) - err(i,j,1) + err(i,j,1)
if ischar(model), model=str2func(model); end
nDataset = size(y, 1); % number of rows, each represents different datasets.
if size(y,2) == 1
    fprintf('X is %d x %d, Y is %d x %d\n',size(x,1),size(x,2),size(y,1),size(y,2));
    warning('It is unlikely you wanted to fit a single data point.  Transposing Y');
    y = transpose(y);
end
if size(x,2) ~= size(y,2)
    fprintf('X is %d x %d, Y is %d x %d\n',size(x,1),size(x,2),size(y,1),size(y,2));
    warning('X and Y have different dimensions.');
end
if size(x, 1) == 1
    x = repmat(x, nDataset, 1);
end
if isa(beta0, 'function_handle'), beta0 = {beta0}; end
if isreal(beta0) && size(beta0, 1) == 1 || length(beta0) == 1
    beta0 = repmat(beta0, nDataset, 1);
end
if exist('mask','var')
    mask = logical(mask);
else 
    mask = [];
end
if isopt(opts,'fine')
    options=statset('TolX',1e-20,'TolFun',1e-20);
else
    options=statset();
end
if isopt(opts, 'robust'), options=statset(options,'Robust','on'); end
if isopt(opts, 'woff')
    ws(1) = warning('query', 'stats:nlinfit:IllConditionedJacobian');
    ws(2) = warning('query', 'stats:nlinfit:IterationLimitExceeded');
    ws2 = ws;
    [ws2.state] = deal('off');
    warning(ws2);
end
for i = 1:nDataset
    if (isopt(opts, 'plfit') || isopt(opts,'plinit')) && ~isopt(opts,'samefig') %  Set up figure 
        figure(500); clf; hold on;
    end 
    if iscell(beta0) % find initial guesses. 
        if isreal(beta0{i})
            beta2 = beta0{i};
        else % assume it's a function
            beta2 = beta0{i}(x(i, :), y(i, :));
        end
    elseif isstruct(beta0)      
        if ~isempty(beta0(i).fn) 
            beta2 = beta0(i).fn(x(i, :), y(i, :), beta0(i).args{:});
            if isfield(beta0, 'vals') % if func and vals, just assign the finite ones as numbers. 
                beta2(isfinite(beta0(i).vals)) = beta0(i).vals(isfinite(beta0(i).vals));
            end
        else % assume you're just using numbers if no function. 
            beta2 = beta0(i).vals;
        end                    
    else
        beta2 = beta0(i, :);
    end
    if i == 1 % check number of params, create parameter list. 
        numParam = length(beta2);        
        if isempty(mask)
            mask = true(1, numParam);
        end
        beta1 = zeros(nDataset, numParam);
    end    
    beta1(i, :) = beta2;        
    plot(x(i, :), y(i, :), '.-'); 
    if isopt(opts, 'plinit') 
        plot(x(i, :), model(beta1(i, :), x(i, :)), 'r--');
    end    
    if ~isopt(opts, 'nofit') % Perform fit. 
      [betaT,rT,JT,COVBT,mseT] = nlinfit(x(i, :), y(i, :), @fitfn, beta1(i, mask),options);
      beta1(i, mask) = betaT;
      res(i,:)=rT(:);
      jac(i,:)=JT(:);
      COVB(i,:)=COVBT(:);
      mse(i,:)=mseT(:);
      try % Find errors. 
        ci = nlparci(betaT,rT,'Jacobian',JT);
        err(i,mask,1:2)=repmat(beta1(i,mask)',1,2)-ci;
      catch
         warning('Could not propagate errors properly')
         err(i,mask,1:2) = nan(sum(mask),2);
      end      
    end  
    if isopt(opts, 'plfit') % Plotted fitted data in red. 
        plot(x(i, :), model(beta1(i, :), x(i, :)), 'r','Linewidth',2);
        xlabel('Frequency (GHz)');
        ylabel('Linearized Power');
        set(gca,'Fontsize',14)
    end    
    if isopt(opts, 'pause') && (i < nDataset), pause; end    
    if isopt(opts,'resid') % Plot residuals. 
      f=gcf; figure(501);
      if ~isopt(opts,'samefig')
        clf;
      else
        hold on;
      end
      plot(x(i,:),y(i,:)-model(beta1(i,:),x(i,:)),'rx-','Color','r','Linewidth',3);
      figure(f);
    end
end   % fit data 
if isopt(opts, 'woff'), warning(ws); end
err = squeeze(err); 

function y = fitfn(beta, x)
beta([find(mask), find(~mask)]) = [beta, beta2(~mask)];
y = model(beta, x);
end
end