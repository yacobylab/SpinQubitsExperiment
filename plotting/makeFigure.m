function f = makeFigure(n)

try
    set(0,'CurrentFigure',n); f = gcf; clf;
catch
    f = figure(n);
end
   
end