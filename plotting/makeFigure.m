function f = makeFigure(n)
% If figure n exists, make it the current figure and clear it, otherwise
% create it. 
% Return figure handle. 
% function f = makeFigure(n)

try
    set(0,'CurrentFigure',n); f = gcf; clf;
catch
    f = figure(n);
end   
end