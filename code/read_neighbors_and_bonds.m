%
% Author: Yosyp Schwab
%
% Function takes x,y,z coordinates and outputs 2 square matrices with
% length equal to number of atoms. First matrix represents all interatomic
% distances between all atoms. Second matrix represents nearest neighbors
% as bonds between atomic pairs.
%
% Inputs:
%        natoms number of atoms <integer>
%        x x-coordinates of atoms <vector>
%        y y-coordinates of atoms <vector>
%        z z-coordinates of atoms <vector>
% Outputs:
%       neighbs interatomic distances between all pairs <sq-matrix>
%       bonds nearest neighbors connected graph <sq-matrix>
%       bond_lengths 1D vector of all bond lengths used for atomic pairs
%
function [neighbs, bonds, bond_lengths] = ...
                                 read_neighbors_and_bonds(natoms, x, y, z)
                             
    bond_lengths = []; 
    for i = 1:natoms
        for j = 1:natoms
            if i ~= j
                neighbs(j,i) = sqrt((x(j)-x(i))^2 + ... 
                                    (y(j)-y(i))^2 + ... 
                                    (z(j)-z(i))^2);
                neighbs(j,i) = round(neighbs(j,i), 4);
            else
                neighbs(j,i) = NaN;
            end
        end
    end

    bonds = zeros(natoms,natoms);
    for i = 1:natoms
        % Alternative approach using:
        % nnindex = find(neighbs(i,:) == min(neighbs(i,:)));

        % Sort all neighbors to find first 3 nearest 
        [out,idx] = sort(neighbs(i,:));

        % In a homogeneus n-ring system only connecting two neighors will
        % form a fully connected graph
        for k = 1:2
            bonds(i,idx(k)) = 1;
            bond_lengths = [bond_lengths out(k)];
        end
        
        % If the system is heterogeneus with 5-6-7-rings, then all three
        % nearest neighbors are sometimes necessary to form a fully
        % connected graph. There is a cutoff (10%) to make sure non-nearest
        % neighbors are not connected.
        if out(3) < 1.3*out(2)
            bonds(i,idx(3)) = 1;
            bond_lengths = [bond_lengths out(3)];
        end
    end
end
