function COP = variancePlot(errorMap, fpMeanCOP,fpVarCOP, stickMeanCOP, stickVarCOP,sphereFit,sphereFitVar)
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
% plot(errorMap(:,1),errorMap(:,2),'xb');
plot(fpMeanCOP(:,1),fpMeanCOP(:,2),'xr');
drawCircle(fpMeanCOP,getRadius(fpMeanCOP, fpVarCOP),'r');
plot(stickMeanCOP(:,1),stickMeanCOP(:,2),'xk');
drawCircle(stickMeanCOP,getRadius(stickMeanCOP, stickVarCOP),'k');
% plot(sphereFit(:,1),sphereFit(:,2),'xg');
% drawCircle(sphereFit,getRadius(sphereFit, sphereFitVar),'g');

% f = zeros(4, 1);
% f(1) = plot(NaN,NaN,'xb');
% f(2) = plot(NaN,NaN,'xr');
% f(3) = plot(NaN,NaN,'xk');
% f(4) = plot(NaN,NaN,'xg');
% legend(f, 'Measured locations','Force plate COP (local coordinates)','Pole tip (rigid body transformation)','Pole tip (sphere fit)');

f = zeros(2, 1);
f(1) = plot(NaN,NaN,'xr');
f(2) = plot(NaN,NaN,'xk');
legend(f,'Force plate COP (local coordinates)','Pole tip (rigid body transformation)');

xlabel('X coordinate(mm)');
ylabel('Y coordinate(mm)');
title('Force plate COP error map')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
yticks(-(h/2):25:(h/2))
xticks(-(w/2):25:(w/2))

    function radius = getRadius(COPmean, COPvar)
       uX = COPmean(:,1);
       uY = COPmean(:,2);
       varX = COPvar(:,1);
       varY = COPvar(:,2);
       %r = sqrt(x^2 + y^2)
       %dr/dx = 0.5*(x^2 + y^2)^(-0.5)*(2x)
       %dr/dy = 0.5*(x^2 + y^2)*(2y)
       %varR = (dr/dx)^2*varX + (dr/dy)^2*varY
       varR = varX.*(0.5 * (uX.^2 + uY.^2).^(-0.5).*(2*uX)).^2 + ...
           varY.*(0.5 * (uX.^2 + uY.^2).^(-0.5).*(2*uY)).^2;
       radius = sqrt(varR);
    end
    function drawCircle(coordinates, radius,color)
        bottomLeft = coordinates - radius;
        topRight = coordinates + radius;
        vector = [bottomLeft radius*2 radius*2];
        for i = 1:length(coordinates)
            rectangle('Position',vector(i,:),'Curvature',[1 1],'EdgeColor',color);
        end
    end
end
