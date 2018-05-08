clear;clc;

path = 'md-runs/rigid-33/stretch_data/';
simfiles = dir(fullfile(path, '*.d'));
nfiles = length(simfiles);

bond_lengths = cell(1);

for i=1:nfiles
    tic
    CNT = importdata(fullfile(path, simfiles(i).name));
    x = CNT(:,3);
    y = CNT(:,4);
    z = CNT(:,5);
    [~,~, bond_lengths{i}] = read_neighbors_and_bonds(length(x), x, y, z);
    fprintf("%d/%d done!\n", i, nfiles);
    toc
end
%%
% natoms = length(x);

% figure;histogram(bond_lengths)
% Ctrs = [0.9:.05:1.805];
% Xcts = hist(bond_lengths{10}, Ctrs);
% bar(Ctrs, Xcts)
% histogram(bond_lengths{1}, .85:.05:1.85)
% histogram(bond_lengths{1}, 20, 'BinLimits', [.8, 1.85])


h = figure;
filename='histogram5945.gif';    
for i = 1:length(bond_lengths)
    histogram(bond_lengths{i}, 1.3:.01:1.75); ylim([0 100]); grid on;
    title(sprintf('Distribution of Bond Lengths, t=%4.4f', i-2));
    xlabel('Bond length [Angstrom]', 'FontWeight', 'bold', 'Color', 'black');
    ylabel('Count [#]', 'FontWeight', 'bold', 'Color', 'black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16); set(gca, 'LineWidth', 2);      
    drawnow;
    
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 

      % Write to the GIF File 
      if i == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
      end         
end