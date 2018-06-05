function COP = arrowPlot(map, fpMeanCOP)
figure('pos',[10 10 1000 1000])
hold on
grid on
w = 500;
h = 300;
if(any(fpMeanCOP(:,2) > h/2))
    h = 600;
end
axis([-(w/2 + 50) (w/2 + 50) -(h/2 + 50) (h/2 + 50)])
rectangle('Position',[-w/2 -h/2 w h], 'EdgeColor','b');
%plot(map(:,1),map(:,2),'xb');

map = map(:,1:2);
dp = fpMeanCOP - map;
quiver(map(:,1),map(:,2),dp(:,1),dp(:,2),0,'Color','r','MaxHeadSize',0.5,'LineWidth', 2.0)

xlabel('X coordinate(mm)');
ylabel('Y coordinate(mm)');
title('COP vector error map')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
yticks(-(h/2):25:(h/2))
xticks(-(w/2):25:(w/2))

end