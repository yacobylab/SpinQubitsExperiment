  function Autoscan()
        
        global tuneData; 
   
        %Run charge scan for autotune    
        tuneData.chrg.run
        
        %Examine the auto fitting to find triple points 
        info1=input('Want to do zoom scan? yes/no:', 's');
        doity1 = strcmp(info1, 'yes');
        doitn1 = strcmp(info1, 'no');
        if doitn1
                tuneData.chrg.ana('man mnslp last'); %manually defining triple points
                tuneData.zoom.run  %Run zoom scan for autotune
              
                info2=input('Accept measpt? yes/no:', 's');
                doity2 = strcmp(info2, 'yes');
                doitn2 = strcmp(info2, 'no');
                
                if doity2
                   tuneData.LoadPos.run;tuneData.LoadTime.run; tuneData.stp.run; tuneData.tl.run 
                end
                
                if doitn2
                    tuneData.zoom.ana('man last') %Manually defining measurement point
                    tuneData.loadPos.run;tuneData.loadTime.run; tuneData.stp.run; tuneData.tl.run 
                end
         end
        
         if doity1
              tuneData.zoom.run  %Run zoom scan for autotune  
              info3=input('Accept measpt? yes/no:', 's');
                doity3 = strcmp(info3, 'yes');
                doitn3 = strcmp(info3, 'no');
                
                if doity3
                   tuneData.loadPos.run; tuneData.loadTime.run; tuneData.stp.run; tuneData.tl.run 
                end
                
                if doitn3
                    tuneData.zoom.ana('man last') %Manually defining measurement point
                    tuneData.loadPos.run;tuneData.loadTime.run; tuneData.stp.run; tuneData.tl.run 
                end
         end
        
end



