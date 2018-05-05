%
% Author: Yosyp Schwab
%
% Function plots <x,y,z> data and optionally shows balls for atoms, sticks
% for bonds, and numbers of atomic indeces.
% Inputs:
%        fig_handle handle to figure for plot()
%        natoms number of atoms <integer>
%        x x-coordinates of atoms <vector>
%        y y-coordinates of atoms <vector>
%        z z-coordinates of atoms <vector>
%        bonds nearest neighbors connected graph <sq-matrix>
%        label_atoms label atoms by number <true/false>
%        ball_size size of balls to displace <float>
%        show_bonds show sticks for bonds <true/false>
%        dskin axes length padding <integer>
% Outputs:
%        none, updated fig_handle only
%
function plot_3d_structure(fig_handle,natoms,x,y,z, ...
                  bonds,label_atoms,ball_size,show_bonds,dskin)

    % Set figure handle
    figure(fig_handle);

    % Set marker size vector
    s = ball_size*ones(1,natoms)';

    hold on;
    grid on;
    view(40,35);
    
%     xlim([min(x)-dskin max(x)+dskin]);
%     ylim([min(y)-dskin max(y)+dskin]);
%     zlim([min(z)-20*(max(x)/max(z)) max(z)+20*(max(x)/max(z))]);
    
    scatter3(x,y,z,s,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1]);

    % Extract atomic pairs from bond matrix, only upper triangular portion
    % otherwise double counting.
    all_atomic_pairs = triu(bonds,1);
    
    for i = 1:natoms
        atom_pair = find(all_atomic_pairs(i,:));
        for k = atom_pair
            if show_bonds
                line([x(i) x(k)], ...
                     [y(i) y(k)], ...
                     [z(i) z(k)], 'Color', 'red', 'LineWidth', 3);
            end
            if label_atoms
                mytxt = text(x(i),y(i),z(i),sprintf('%d',i));         
                mytxt.Color = 'blue';
                mytxt.FontSize = 16;
            end
        end
    end
    
end
