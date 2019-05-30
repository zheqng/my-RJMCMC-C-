
delete canadianweather_temp.dat
fid = fopen('CanadianWeather_temperature.dat','w')
for ii = 1:17
    for jj = 1:37
        fprintf(fid,'%f ',T1{1,ii}(jj));
    end
    fprintf(fid,'\n');
     fprintf(fid,'\n');
  for jj = 1:37
        fprintf(fid,'%f ',Y1{1,ii}(jj));
    end
    fprintf(fid,'\n');
     fprintf(fid,'\n');
end
fclose(fid);
