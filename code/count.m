%
% To make CNT's see:
% http://turin.nss.udel.edu/research/tubegenonline.html
% To post-process the files see:
% convert_cnt_coords.sh
%
clear;clc;

% CNT = dlmread('../input_files/cnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/longcnt.xyz', ' ');
% CNT = dlmread('../input_files/midcnt.xyz', ' ');
% CNT = dlmread('../input_files/cappedcnt.xyz', ' ');
% CNT = dlmread('../input_files/gencnt.xyz', ' ');
% CNT = dlmread('../input_files/c60.xyz');
% CNT = dlmread('../input_files/cnt-33-112.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cnt-33-1110-defect.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cnt-33-1110.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cnt-1010-1110.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cnt-170-1110.xyz', ' ', 2, 1);
% 
% x = CNT(:,1);
% y = CNT(:,2);
% z = CNT(:,3);
% natoms = length(x);

% CNT = importdata('md-runs/data/time00100.d');
% % CNT = importdata('md-runs/equilibration_data/time00001.d');
CNT = importdata('md-runs/rigid-170/stretch_data/time00047.d');
x = CNT(:,3);
y = CNT(:,4);
z = CNT(:,5);
natoms = length(x);

[neighbs, bonds, bond_lengths] = read_neighbors_and_bonds(natoms, x, y, z);
% figure;histogram(bond_lengths)

pos_fig     = figure('Position', [1000 100 500 1200]);
ball_size   = 250;
label_atoms = false;
show_bonds  = true;
dskin       = 0;
plot_3d_structure(pos_fig,natoms,x,y,z, ...
                  bonds,label_atoms,ball_size,show_bonds,dskin);

% bonds = (bonds + bonds')/2;              
atom_graph = graph(bonds);

%     fig_graph = figure;
%         h=plot(atom_graph);
%         title('Undirected Graph of System');
%         xt = get(gca, 'XTick'); set(gca, 'FontSize', 16); set(gca, 'LineWidth', 2);        
% highlight(h,[23,24],'NodeColor','g','EdgeColor','g','MarkerSize',15,'LineWidth',4)
% highlight(h,[35,24],'NodeColor','g','EdgeColor','g','MarkerSize',15,'LineWidth',4)
% highlight(h,[25,24],'NodeColor','g','EdgeColor','g','MarkerSize',15,'LineWidth',4)
% 
% highlight(h,[35,36],'NodeColor','r','EdgeColor','r','MarkerSize',15,'LineWidth',4)
% highlight(h,[35,34],'NodeColor','r','EdgeColor','r','MarkerSize',15,'LineWidth',4)
% highlight(h,[25,26],'NodeColor','r','EdgeColor','r','MarkerSize',15,'LineWidth',4)
% highlight(h,[25,14],'NodeColor','r','EdgeColor','r','MarkerSize',15,'LineWidth',4)
% highlight(h,[23,12],'NodeColor','r','EdgeColor','r','MarkerSize',15,'LineWidth',4)
% highlight(h,[23,22],'NodeColor','r','EdgeColor','r','MarkerSize',15,'LineWidth',4)
% 
% highlight(h,[34,33],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[22,33],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[22,21],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[12,13],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[12,11],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[14,13],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[14,15],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[26,27],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[26,37],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)
% highlight(h,[36,37],'NodeColor','k','EdgeColor','k','MarkerSize',15,'LineWidth',4)


figure(pos_fig);
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated34asdasdfsadfadssadf22.gif';
ax = gca;
ax.CameraTargetMode = 'manual';
ax.CameraTarget = [ mean(x) mean(y) mean(z) ];
ax.CameraViewAngleMode = 'manual';
ax.CameraViewAngle = 100;
theta = 0;
for theta=0:-.02:-2*pi
    ax.CameraPositionMode = 'manual';
    ax.CameraPosition = [ mean(x)+5*cos(theta) mean(y)+5*sin(theta) mean(z)+0 ];    
    ax.XLim = [min(z)-5 max(z)+5];
    ax.YLim = [min(z)-5 max(z)+5];
    ax.ZLim = [min(z)-5 max(z)+5];

    ax.Projection = 'orthographic';
    grid off;
    ax.Color = 'none';
    ax.Visible = 'off';
    drawnow;
      % Capture the plot as an image 
      frame = getframe(pos_fig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 

      % Write to the GIF File 
      if theta == 0
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
      end     
%     pause(.04);
end
%%
% figure(pos_fig);
% ax = gca;
% ax.CameraPositionMode = 'manual';
% ax.CameraPosition = [ mean(x) mean(y) mean(z)+0 ];
% ax.CameraViewAngleMode = 'manual';
% ax.CameraViewAngle = 120;
% for theta=0:-.1:-2*pi
%     ax.CameraTargetMode = 'manual';
%     ax.CameraTarget = [ 1*cos(theta)-mean(x) 1*sin(theta)+mean(y) mean(z)+0 ];
%     ax.XLim = [min(z) max(z)];
%     ax.YLim = [min(z) max(z)];
%     ax.ZLim = [min(z) max(z)];
%     ax.Projection = 'orthographic';
%     drawnow;
%     pause(.04);
% end



%%
% [fig_bonds,fig_neighbs,fig_graph] = plot_info(bonds,neighbs, atom_graph);
%%
% first_node = 2;
first_node = 24;
% first_node = 24;

[first_neighbors, ...
 second_neighbors, ...
 third_neighbors] = find_extended_neighbors(atom_graph, first_node);

% Find all connected 4-segments stemming from first node
% Store all found segments in rings matrix
segments  = zeros(3,natoms);
five_seg  = zeros(3,natoms);
six_seg   = zeros(3,natoms);
seven_seg = zeros(3,natoms);
counter   = 1; % Count of total found 4-segments
for i=1:length(first_neighbors)
    for j=1:2 % second_neighbor row counter
        for k=1:2 % third_neighbor column counter
            
            % Each second_neighbor row has two atoms, for a total of 4
            % third_neighbor atoms. The third_neighbor row counter should
            % update every other k-loop.
            if j == 1
                n3k = 2*i-1; % third_neighbor row counter
            else
                n3k = 2*i; % third_neighbor row counter
            end
            
            if (isnan(third_neighbors(n3k,k)) == false) && (isnan(second_neighbors(i,j)) == false)
%                 fprintf("Loop:%d-%d-%d-%d: %d -> %d -> %d -> %d\n", ... 
%                        counter, i,j,k, first_node, ...
%                        first_neighbors(i), second_neighbors(i,j), third_neighbors(n3k,k));
                fprintf("Connected Segment: %d -> %d -> %d -> %d\n", ... 
                        first_node, ...
                       first_neighbors(i), second_neighbors(i,j), third_neighbors(n3k,k));
                seven_seg(counter, third_neighbors(n3k,k)) = 1;
                six_seg(counter,[first_node, ...
                                 third_neighbors(n3k,k)]) = 1;
                five_seg(counter,[first_node, ...
                                  second_neighbors(i,j), ...
                                  third_neighbors(n3k,k)]) = 1;                             
                segments(counter,[first_node, ...
                                  first_neighbors(i), ...
                                  second_neighbors(i,j), ...
                                  third_neighbors(n3k,k)]) = 1;
                counter = counter+1;    
            end
        end
    end
end
%%
twopairs = nchoosek(1:counter-1, 2);
five_rings  = [];
six_rings   = [];
seven_rings = [];
for i=1:length(twopairs)

    % 6-rings
    % Find number of common nodes between all pairs of six_segments which
    % only have at two nodes
    six_pairs =  sum(bitand( ...
                             six_seg(twopairs(i,1),:), ...
                             six_seg(twopairs(i,2),:) ) );
    % six_seg's with two nodes identical first and third neighbor nodes
    % form a ring. Store the full node-member list of the ring.
    if six_pairs == 2
         six_rings(end+1,:) = find(bitor( ...
                                          segments(twopairs(i,1),:), ...
                                          segments(twopairs(i,2),:) ));
    end
    
%     if bonds( ) == 1
%          six_rings(end+1,:) = find(bitor( ...
%                                           segments(twopairs(i,1),:), ...
%                                           segments(twopairs(i,2),:) ));
%     end
    
    % 5-rings
    five_pairs = sum(bitand( ...
                             five_seg(twopairs(i,1),:), ...
                             five_seg(twopairs(i,2),:) ) );                            
    if five_pairs == 3
        five_rings(end+1,:) = find(bitor( ...
                                          segments(twopairs(i,1),:), ...
                                          segments(twopairs(i,2),:) ));
    end
    
    % 7-rings
    s1 = find(seven_seg(twopairs(i,1),:))
    s2 = find(seven_seg(twopairs(i,2),:))
    if bonds(s1,s2) == 1
            fprintf('match!\n');
            segments(twopairs(i,1),:)
            segments(twopairs(i,2),:)
            bitor(segments(twopairs(i,1),:),segments(twopairs(i,2),:))
            find(bitor(segments(twopairs(i,1),:),segments(twopairs(i,2),:)))
          seven_rings(end+1,:) = find(bitor( ...
                                             segments(twopairs(i,1),:), ...
                                             segments(twopairs(i,2),:) ));
    end

end

figure(pos_fig);
[n6,~] = size(six_rings);
for i=1:n6
    for j=1:6
        scatter3(x(six_rings(i,j)),y(six_rings(i,j)),z(six_rings(i,j)),250,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .75 .75]);
    end
end

[n5,~] = size(five_rings);
for i=1:n5
    for j=1:5
        scatter3(x(five_rings(i,j)),y(five_rings(i,j)),z(five_rings(i,j)),250,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63]);
    end
end

fprintf('7-rings:');
disp(six_rings);
fprintf('6-rings:');
disp(six_rings);
fprintf('5-rings:');
disp(five_rings);

