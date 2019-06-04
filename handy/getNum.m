function num = getNum(str)
num = str2double(regexp(str,'\d+','match'));
end
