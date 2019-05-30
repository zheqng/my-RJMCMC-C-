function plot_curve(T,Y,group)
train_size=numel(T);
figure
hold on
if nargin==3
color={'r','g','b','m','c','k'};
for i=1:train_size
    plot(T{i},Y{i},color{group(i)})
end
else
    for i=1:train_size
        plot(T{i},Y{i})
    end
end
hold off