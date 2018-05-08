%
% To make CNT's see:
% http://turin.nss.udel.edu/research/tubegenonline.html
% To post-process the files see:
% convert_cnt_coords.sh
%
clear;clc;


% CNT = dlmread('../input_files/c60.xyz', ' ');
% CNT = dlmread('../input_files/cappedcnt.xyz', ' ');
% CNT = dlmread('../input_files/longcnt.xyz', ' ');
% CNT = dlmread('../input_files/midcnt.xyz', ' ');
% CNT = dlmread('../input_files/cnt-33-112.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cnt-33-1110-defect.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cnt-33-1110.xyz', ' ', 2, 1);
% x = CNT(:,1);
% y = CNT(:,2);
% z = CNT(:,3);
% natoms = length(x);

% CNT = importdata('md-runs/rigid-33/stretch_data/time00065.d');
CNT = importdata('md-runs/rigid-170/stretch_data/time00047.d');
x = CNT(:,3);
y = CNT(:,4);
z = CNT(:,5);
natoms = length(x);

[neighbs, bonds, bond_lengths] = read_neighbors_and_bonds(natoms, x, y, z);
figure;histogram(bond_lengths)

pos_fig     = figure('Position', [1000 100 500 1200]);
ball_size   = 250;
label_atoms = false;
show_bonds  = true;
dskin       = 0;
plot_3d_structure(pos_fig,natoms,x,y,z, ...
                  bonds,label_atoms,ball_size,show_bonds,dskin);

bonds = (bonds + bonds')/2;              
atom_graph = graph(bonds);

% [fig_bonds,fig_neighbs,fig_graph] = plot_info(bonds,neighbs, atom_graph);


% first_node = 2;
first_node = 1;
% first_node = 24;
five_rings  = [];
six_rings   = [];
seven_rings = [];

for first_node = 1:length(x)
    tic
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
                    %fprintf("Loop:%d-%d-%d-%d: %d -> %d -> %d -> %d\n", ... 
                    %        counter, i,j,k, first_node, ...
                    %        first_neighbors(i), second_neighbors(i,j), third_neighbors(n3k,k));
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

    if length(1:counter-1) > 2
        twopairs = nchoosek(1:counter-1, 2);
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

            % 5-rings
            five_pairs = sum(bitand( ...
                                     five_seg(twopairs(i,1),:), ...
                                     five_seg(twopairs(i,2),:) ) );                            
            if five_pairs == 3
                five_rings(end+1,:) = find(bitor( ...
                                                  segments(twopairs(i,1),:), ...
                                                  segments(twopairs(i,2),:) ));
            end
        end
    end
    toc
end

unique_six_rings = unique(six_rings, 'rows');
unique_five_rings = unique(five_rings, 'rows');

fprintf('Found %d unique 6-rings, they are:\n', length(unique_six_rings(:,1)));
disp(unique_six_rings);
fprintf('Found %d unique 5-rings, they are:\n', length(unique_five_rings(:,1)));
disp(unique_five_rings);


%%
figure(pos_fig);
[n6,~] = size(unique_six_rings);
for i=1:n6
    for j=1:6
        scatter3(x(unique_six_rings(i,j)),y(unique_six_rings(i,j)),z(unique_six_rings(i,j)),1,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .75 .75]);
    end
end

[n5,~] = size(unique_five_rings);
for i=1:n5
    for j=1:5
        scatter3(x(unique_five_rings(i,j)),y(unique_five_rings(i,j)),z(unique_five_rings(i,j)),250,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63]);
    end
end

