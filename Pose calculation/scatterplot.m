COP_error = []
angles = [];
force = [];
combined = [];

[file,path] = uigetfile('*.mat','Select input file(s)','MultiSelect','on');
if(~iscell(file))
    file = {file};
end

recording = {};
rec = {};
for i = 1:length(file)
    combined = [combined load(strcat(path,file{i}))];
end

for i = 1:length(combined)
    force = [force combined(i).errorData.force]
    COP_error = [COP_error combined(i).errorData.COP_error];
    angles = [angles combined(i).errorData.angles];
end

%% 
norm_force = vecnorm(force);
norm_cop = vecnorm(COP_error);
idx = find(vecnorm(COP_error) < 100);
figure('Name','scatter_force')
hold on
scatter(norm_force(idx),norm_cop(idx), 'rx')
xlabel('Force [N]')
ylabel('COP Error [mm]')
title('Scatterplot of Force vs COP error')
[b1, bint1, r1, rint1, stats1] = regress(norm_cop(idx),[ones(size(norm_force(idx))) norm_force(idx)]);
plot([0 max(norm_force(idx))],[b1(1) max(norm_force(idx))*b1(2) + b1(1)],'b')
lm = fitlm(norm_force(idx),norm_cop(idx),'linear')
legend('COP error [mm]', sprintf('Fit line: y = %6.1fx + %6.3f',b1(2), b1(1)))

savefig(strcat(path,'/scatter_force'));

figure('Name','scatter_inv_force')
scatter(1./norm_force(idx),norm_cop(idx), 'rx')
xlabel('Force^{-1} [1/N]')
ylabel('COP Error [mm]')
title('Scatterplot of 1/Force vs COP error')
hold on
[b, bint, r, rint, stats] = regress(norm_cop(idx),[ones(size(1./norm_force(idx))) 1./norm_force(idx)]);
plot([0 max(1./norm_force(idx))],[b(1) max(1./norm_force(idx))*b(2) + b(1)],'b')
lm2 = fitlm(1./norm_force(idx),norm_cop(idx),'linear')
legend('COP error [mm]', sprintf('Fit line: y = %6.1fx - %6.3f',b(2), abs(b(1))))
savefig(strcat(path,'/scatter_inv_force'));
%% 
x1 = linspace(0,max(1./norm_force(idx)));
[p, S] = polyfit(1./norm_force(idx),norm_cop(idx),2);
y1 = polyval(p,x1);

figure 
hold on
scatter(1./norm_force(idx),norm_cop(idx), 'rx')
plot(x1, y1,'b')
lm3 = fitlm(1./norm_force(idx),norm_cop(idx),'quadratic')
save(strcat(path,'/scatterData.mat'))
xlabel('Force^{-1} [1/N]')
ylabel('COP Error [mm]')
title('Scatterplot of 1/Force vs COP error')
legend('COP error [mm]', sprintf('Fit line: y = %3.2fx^{2} - %3.3fx + %6.3f',p(1), abs(p(2)),  p(3)))
savefig(strcat(path,'/scatter_inv_force_poly'));
