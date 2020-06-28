load = zeros(1000,96);
time = repmat(1:96,1000,1);
for ii=1:1000
   index =  ((ii-1)*96+1):((ii-1)*96+96);
    load(ii,:) = result_out(index);
    plot(time(ii,:),load(ii,:));
    hold on
end
hold off


figure
for ii=1:100
    plot(time(ii,:),load(ii,:)/1000);
    hold on
end
hold off

delete('../demo/electric_load_2008.dat');
fid = fopen('../demo/electric_load_2008.dat','w')
for ii = 1:100
    for jj = 1:96
        fprintf(fid,'%f ',time(ii,jj));
    end
    fprintf(fid,'\n');
    % fprintf(fid,'\n');
  for jj = 1:96
        fprintf(fid,'%f ',load(ii,jj)/1000);
    end
    fprintf(fid,'\n');
   %  fprintf(fid,'\n');
end
fclose(fid);


delete('../demo/testdata.dat');
fid = fopen('../demo/testdata.dat','w')

   
for ii = 101:114
  for jj = 1:96
        fprintf(fid,'%f ',load(ii,jj)/1000);
    end
end
    fprintf(fid,'\n');
fclose(fid);
