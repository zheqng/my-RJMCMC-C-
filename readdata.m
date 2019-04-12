for(m = 1:9)
dlmwrite('traindata.dat',(traindata{m}(:,1))','-append','delimiter',' ','roffset',1,'precision',4)
dlmwrite('traindata.dat',(traindata{m}(:,2))','-append','delimiter',' ','roffset',1,'precision',4)
end

for(m = 1:9)
dlmwrite('testdata.dat',(trainrestdata{m}(:,1))','-append','delimiter',' ','roffset',1,'precision',4)
dlmwrite('testdata.dat',(trainrestdata{m}(:,2))','-append','delimiter',' ','roffset',1,'precision',4)
end