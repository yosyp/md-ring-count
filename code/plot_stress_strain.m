clear;clc;

path = 'md-runs/rigid-170/stretch_data/';
simfiles = dir(fullfile(path, '*.d')); % only first 65
nfiles = length(simfiles);
strain = [];
stress = [];

for i=1:nfiles
    tic
    CNT = importdata(fullfile(path, simfiles(i).name));
    strain1(i) = max(CNT(CNT(:,2)~=3,5)) - min(CNT(CNT(:,2)~=3,5));
    stress1(i) = (sum(CNT(CNT(:,2)~=3,6)) + sum(CNT(CNT(:,2)~=3,7)) + sum(CNT(CNT(:,2)~=3,8)))/3;
    fprintf("1: %d/%d done!\n", i, nfiles);
    toc
end

path = 'md-runs/rigid-1010/stretch_data/';
simfiles = dir(fullfile(path, '*.d')); % only first 65
nfiles = length(simfiles);
strain = [];
stress = [];

for i=1:nfiles
    tic
    CNT = importdata(fullfile(path, simfiles(i).name));
    strain2(i) = max(CNT(CNT(:,2)~=3,5)) - min(CNT(CNT(:,2)~=3,5));
    stress2(i) = (sum(CNT(CNT(:,2)~=3,6)) + sum(CNT(CNT(:,2)~=3,7)) + sum(CNT(CNT(:,2)~=3,8)))/3;
    fprintf("2: %d/%d done!\n", i, nfiles);
    toc
end

path = 'md-runs/rigid-33/stretch_data/';
simfiles = dir(fullfile(path, '*.d')); % only first 65
nfiles = length(simfiles);
strain = [];
stress = [];

for i=1:nfiles
    tic
    CNT = importdata(fullfile(path, simfiles(i).name));
    strain3(i) = max(CNT(CNT(:,2)~=3,5)) - min(CNT(CNT(:,2)~=3,5));
    stress3(i) = (sum(CNT(CNT(:,2)~=3,6)) + sum(CNT(CNT(:,2)~=3,7)) + sum(CNT(CNT(:,2)~=3,8)))/3;
    fprintf("3: %d/%d done!\n", i, nfiles);
    toc
end

%%
md2newton = 1.660539E-13;
md2pascal = 1.660539E7;
% md2pascal = 1;

k = 47;
normstrain1 = ((strain1(1:k)-strain1(1))/strain1(1));
normstress1 = (-stress1(1:k)+stress1(1)).*md2pascal./1e12;

k = 54;
normstrain2 = ((strain2(1:k)-strain2(1))/strain2(1));
normstress2 = (-stress2(1:k)+stress2(1)).*md2pascal./1e12;

k = 63;
normstrain3 = ((strain3(1:k)-strain3(1))/strain3(1));
normstress3 = (-stress3(1:k)+stress3(1)).*md2pascal./1e12;

plot(normstrain1, normstress1, ...
     normstrain2, normstress2, ...
     normstrain3, normstress3, 'LineWidth', 5); grid on;
title('Stress vs Strain @ Different Chiralities');
legend('(17,0)','(10,10)','(3,3)');
xlabel('Strain','FontWeight', 'bold', 'Color', 'black');
ylabel('Stress [GPa]','FontWeight', 'bold', 'Color', 'black');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16); set(gca, 'LineWidth', 2);        
    
P = polyfit(normstrain3,normstress3,1);
% yfit = P(1)*normstrain;
disp(P(1));
% plot(normstrain,normstress); hold on;
% xlabel('Strain ');
% ylabel('Stress [GPa]');
% plot(normstrain,yfit, 'r-.');
% disp(P(1))
