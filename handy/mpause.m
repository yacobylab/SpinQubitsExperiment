function mpause(x)
% More accurate for short pauses than Matlab

if(x > 0.01)
pause(x);
else
    tic;while(toc < x) ; end; 
end
