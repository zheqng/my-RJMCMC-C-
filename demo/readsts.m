clear all;clc;close all;
fid = fopen('CanadianWeather_temperatureRes.sts','r');

n = 10000;
M = 3;
X = zeros(1,10);Y = zeros(n,10);Class = [];
W = [];V0 = []; Sigmav2 = [];Pi = [];k=8;neg_count=0;
for i=1:n
    iter(i)=fscanf(fid,'%dth iter\n',[1 1]);
    K(i) = fscanf(fid,'%d\n',[1 1]);
%     fscanf(fid,'%f',[1 1]);
%  fscanf(fid,'theta->pi is:');
 PI{i} = fscanf(fid,'%f',[1 K(i)]);
    w{i} = fscanf(fid,'%f',[1 K(i)]);
   
    v0{i}=fscanf(fid,'%f',[1 K(i)]);
%     a0(i)=fscanf(fid,'%f',[1 1]);
%     a1(i) = fscanf(fid,'%f',[1 K(i)]);
    sigmav2{i} = fscanf(fid,'%f',[1 K(i)]);
     if sum(v0{i}<0)>0
        neg_count = neg_count+1;
    end
    Z(i,:)=fscanf(fid,'%d',[1 20]);
%     X(K(i):end) = X(K(i):end)+1;
  X(K(i)) = X(K(i))+1;
    Y(i,:) = X/i;
    
    if K(i) ==k
        Pi = [Pi;PI{i}];
        W = [W;w{i}];
        V0 = [V0;v0{i}];
        Sigmav2 = [Sigmav2;sigmav2{i}];
        Class = [Class ;Z(i,:)];
    end
end
%____________________hist K_______________________%
x=[2 3 4  5 6];
hist(K/n,x)

figure
plot(Y)
figure
for i =1:6
subplot(2,3,i)
plot(Y(:,i));hold on
end

figure
bar(Y(end,:))
%___________________kernel smoothing density,posterior density of estimates_______________%

%_____________________pi_______________________%
figure
% x=-2000:-1000;

x =0:0.001:1;
for i=1:k
    pd = fitdist(Pi(:,i),'kernel','Kernel','normal');
%  pd = fitdist(-2*L{i}','gamma');
% pd = chi2pdf(-2*);
  
    y = pdf(pd,x);
    subplot(4,k,i);
    plot(x,y);
    title('pi')
end
%____________________________w_______________%
% figure
% x=-2000:-1000;

x = 0:0.001:2;
for i=1:k
    pd = fitdist(W(:,i),'kernel','Kernel','normal');
%  pd = fitdist(-2*L{i}','gamma');
% pd = chi2pdf(-2*);
  
    y = pdf(pd,x);
    subplot(4,k,k+i);
    plot(x,y);
    title('w')
end
% x = 0:0.001:8;
%  pd = fitdist(W(:,3),'kernel','Kernel','normal');
% %  pd = fitdist(-2*L{i}','gamma');
% % pd = chi2pdf(-2*);
%   
%     y = pdf(pd,x);
%     subplot(4,k,k+3);
%     plot(x,y);
%     title('w')
%_____________________v0_______________________%
% figure
% x=-2000:-1000;

% x =0:0.001:2;
% for i=1:k
%     pd = fitdist(V0(:,i),'kernel','Kernel','normal');
% %  pd = fitdist(-2*L{i}','gamma');
% % pd = chi2pdf(-2*);
%   
%     y = pdf(pd,x);
%     subplot(4,k,2*k+i);
%     plot(x,y);
%     title('v0')
% end
% x =0:0.001:4;
%  pd = fitdist(V0(:,2),'kernel','Kernel','normal');
% %  pd = fitdist(-2*L{i}','gamma');
% % pd = chi2pdf(-2*);
%   
%     y = pdf(pd,x);
%     subplot(4,k,2*k+2);
%     plot(x,y);
%     x =0:0.001:8;
%  pd = fitdist(V0(:,1),'kernel','Kernel','normal');
% %  pd = fitdist(-2*L{i}','gamma');
% % pd = chi2pdf(-2*);
%   
%     y = pdf(pd,x);
%     subplot(4,k,2*k+1);
%     plot(x,y);
% % x = 0:0.001:25;
% % pd = fitdist(V0(:,2),'kernel','Kernel','normal');
% % y = pdf(pd,x);
% % subplot(4,k,2*k+2);
% %   plot(x,y);
% % title('v0');
% %_____________________sigmav2_______________________%
% % figure
% x =0:0.00001:0.001;
% for i=1:k
%     pd = fitdist(Sigmav2(:,i),'kernel','Kernel','normal');
% %  pd = fitdist(-2*L{i}','gamma');
% % pd = chi2pdf(-2*);
%   
%     y = pdf(pd,x);
%     subplot(4,k,3*k+i);
%     plot(x,y);
% %     ylim([0,35]);
%     title('sigmav2')
% end
% x =0:0.00001:0.005;
%   pd = fitdist(Sigmav2(:,2),'kernel','Kernel','normal');
% %  pd = fitdist(-2*L{i}','gamma');
% % pd = chi2pdf(-2*);
%   
%     y = pdf(pd,x);
%     subplot(4,k,3*k+2);
%     plot(x,y);
% %     ylim([0,35]);
%     title('sigmav2')
% % x = -0.001:0.000001:0.001;
% % pd = fitdist(Sigmav2(:,2),'kernel','Kernel','normal');
% % %  pd = fitdist(-2*L{i}','gamma');
% % % pd = chi2pdf(-2*);
% %   
% %     y = pdf(pd,x);
% %     subplot(4,k,3*k+2);
% %     plot(x,y,'.');
% % %     ylim([0,3.5]);
% %     title('sigmav2')
% %     figure
% %     hist(Sigmav2(:,2))

for i = 1:k
   PImode(i) =  mode(Pi(:,i));
   Wmode(i) = mode(W(:,i));
   V0mode(i) = mode(V0(:,i));
   Sigmav2mode(i) = mode(Sigmav2(:,i));
end
PImode
Wmode
V0mode
Sigmav2mode
neg_count



figure
label = mode(Class,1);
label = label +1;
% colors = {'r'  'g' 'b' 'm' 'cyan' 'yellow' [1 0 0]};
% ∫Ï¬Ã¿∂«‡—Û∫Ïª∆≥»
c = [1 0 0;0 1 0; 0 0 1;0 1 1;1 0 1;1 1 0; 1 0.5 0;0.5 0 0];
lines = {'-','-','-.','--',':',':'};
load paper.mat
% plot(agefine, hgtfmat2)
% fid = fopen('girlsgrowthcurves.dat','r');
% A = fscanf(fid,'%f ',[ 8 12*30]);
% data  = cell(30,1);
% gradient = cell(30,1);
[Nm,Curve_num] = size(hgtfmat2);
% Curve_num = size(gradient,1);
% Nm = size(gradient{1},2);
figure
for m = 1:Curve_num
    plot(agefine,hgtfmat2(:,m),'Color',c(label(m),:));hold on
end
hold off
% % real cluster
% atlindex = [1,2,3,4,5,6,7,8,9,10,11,13,14,16];
% pacindex = [25,26,27,28,29];
% conindex = [12,15,17,18,19,20,21,22,23,24,30,31,35];
% artindex = [32,33,34];
% real_label = zeros(1,35);
% real_label(atlindex)=1;
% real_label(pacindex)=2;
% real_label(conindex)=3;
% real_label(artindex)=4;
% figure
% for m = 1:Curve_num
%     plot(monthtime(1:end-1),gradient(:,m),'Color',colors{real_label(m)},'LineStyle',lines{real_label(m)});hold on
% end
% hold off

