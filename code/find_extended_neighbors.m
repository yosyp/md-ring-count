%
% Author: Yosyp Schwab
%
% Function plots <x,y,z> data and optionally shows balls for atoms, sticks
% for bonds, and numbers of atomic indeces.
% Inputs:
%        atom_graph undirected graph with nodes = coordinates <graph>
%        first_node atom index where to begin searching <integer>
% Outputs:
%        first_neighbors nearest-neighbors <3x1 vector>
%        second_neighbors nearest-nearest-neighbors <3x2 matrix>
%        third_neighbors nearest-nearest-nearest-neighbors <6x2 matrix>
%
function [first_neighbors, second_neighbors, third_neighbors] = ...
                            find_extended_neighbors(atom_graph,first_node)

    % Using the first_node, first_neighbors will be found consisting of
    % three nodes <3x1 vector>.
    %
    % For each first_neighbor node, two nearest neighbor nodes will be
    % found (not backlinking to original node) <3x2 matrix> with each row
    % corresponding to the originating first_neighbor.
    %
    % For each second_neighbor, two nearest neighbors nodes will be found
    % (not backlinking to original node) <6x2 matrxi> with each row
    % correspoding to the originating second_neighbor.
  
    % Edge-nodes sometimes do not have 2 neighbor nodes, skips and
    % matrix-padding are required (NaN used for missing node).
    
    % The structure/relation/alignmemnt of the three datasets are important
    % because that is how segments/loops are found.
    
    second_neighbors = NaN(3,2);
    third_neighbors  = NaN(6,2);
    first_neighbors  = neighbors(atom_graph,first_node);
    
    k=1;
    for i=1:length(first_neighbors)
        
        % Find second neighbors, remove backlink to the "calling node", and
        % pad the row if necessary (if less than 2 neighbor nodes found).
        % Store processed row in second_neighbors.
        nb2 = neighbors(atom_graph,first_neighbors(i));
        nb2 = nb2(nb2 ~= first_node);
        if length(nb2) < 2
            nb2 = [nb2; NaN];
        end    
        
        second_neighbors(i,:) = nb2;
        
        for j=1:length(second_neighbors(i,:))
            if isnan(second_neighbors(i,j)) == true
                third_neighbors(k,:) = [NaN NaN];
            else
                nb3 = neighbors(atom_graph, second_neighbors(i,j));
                nb3 = nb3(nb3 ~= first_neighbors(i));  
                if length(nb3) < 2
                    nb3 = [nb3; NaN];
                end
                third_neighbors(k,:) = nb3;
            end
            k=k+1;
        end

    end


end
