function make2grp(pg1,pg2,name)
% Simplify making 2 qubit groups. Only works for typical types. 
pg.ctrl='grp loop pack';
npulses = size(pg1.varpar,1); 
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pg1.chan, pg2.chan];
pg.matrix=eye(length(pg.chan));

pg.pulses.groups={pg1.name,pg2.name};
pg.name=name;
plsdefgrp(pg);
end