

clear all;clc;close all;
fid = fopen('CanadianWeather_temperature.sts','r');
%  a=fscanf(fid,'%g%g',[9 inf]);

% n=20000;
% X=[];Y=[];
% for i=1:n
%    a = fscanf(fid,'%d %d',[1 2]);
%    logl(i) = fscanf(fid,'%f',[1 1]);
%    b = fscanf(fid,'%f\n',[1 1]);
%    if mod(i,20)==0
% %        plot(i,logl(i),'.');hold on;
% X = [X i];Y = [Y logl(i)];
%    end
% end
% plot(X,Y,'.')


n=10000;X=[];Y=[];Z=[];L=cell(10,1);split_count = 0;merge_count = 0;
for i=1:n
   a = fscanf(fid,'%d',[1 1]);
   K(i) = fscanf(fid,'%f',[1 1]);
   logl(i) = fscanf(fid,'%f',[1 1]);
   move = fscanf(fid,'%s',[1 1]);
   accept = fscanf(fid,'%s',[1 1]);
   Paccept(i) = fscanf(fid,'%f',[1 1]);
   c = fscanf(fid,'%s',[1 2]);
%    L{K(i)} = [L{K(i)} logl(i)];
%    b = fscanf(fid,' %f\n',[1 1]);
%    if mod(i,200)==0
%        plot(i,logl(i),'.');hold on;
       L{K(i)} = [L{K(i)} logl(i)];
       if (strcmp(move , 'split') && strcmp(accept , 'accept'))
           split_count = split_count +1;
       end
       if strcmp(move ,'merge') && strcmp(accept,'accept')
           merge_count = merge_count +1;
       end
%       X = [X i];Y = [Y logl(i)];Z=[Z K(i)];
%    end
end
figure
plot(logl,'.');
% plot(K);hold off
% figure
% AX = plotyy(X,Y,X,Z);
% ylim(AX(1),[550,1000]);
% ylim(AX(2),[0 7])
figure
[AX, H1, H2] = plotyy(1:n,K,1:n,Paccept);
ylim(AX(1),[1,9]);
ylim(AX(2),[0,1]);
set(H1,'LineStyle','--')
set(H2,'LineStyle','--')
figure
for i =2:8
    subplot(3,3,i)
    plot(-2*L{i},'.')
end
title("L_i")
%___________________kernel smoothing density_______________%
figure
x=-14500:-12500;
colors = {'g' 'r' 'g' 'b' 'g' 'b' 'm' 'cyan'};
lines = {'-','-','-.','-.','--','--',':','-'};
i=5;
pd = fitdist(-2*L{i}','kernel','Kernel','normal');
%  pd = fitdist(-2*L{i}','gamma');
% pd = chi2pdf(-2*);
    y = pdf(pd,x);
    plot(x,y,'Color',colors{i},'LineStyle',lines{i});hold on;
for i=1:6
    pd = fitdist(-2*L{i}','kernel','Kernel','normal');
%  pd = fitdist(-2*L{i}','gamma');
% pd = chi2pdf(-2*);
    y = pdf(pd,x);
    plot(x,y,'Color',colors{i},'LineStyle',lines{i});hold on;
end
title("fit density")
hold off;
split_count/n
merge_count/n

accept= cell(1,10);
for ii = 1:n
    if(Paccept(ii)>0.5)
    accept{K(ii)} = [accept{K(ii)} Paccept(ii)];
    end
end
figure
for ii = 2:6
    subplot(2,3,ii)
    plot(accept{ii},'.')
end
title("acceptance probability")
