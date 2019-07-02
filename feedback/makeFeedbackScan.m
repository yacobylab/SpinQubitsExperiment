function fbScan = makeFeedbackScan(grp,nloop,datachan)

fbScan=fConfSeq(grp,struct('nloop',nloop,'nrep',1,'opts','raw','datachan',datachan));
fbScan.configch=[];
fbScan.figure=1111;
fbScan.loops(1).setchan='count2'; % Use count2 in case count being set in main scan.
fbScan.disp=fbScan.disp([fbScan.disp.dim]==1); % Only plot 1d data.
fbScan.consts(1) = []; % remove clock setting.
fbScan.xv=fbScan.data.pulsegroups(1).varpar(:,1)';
            
end