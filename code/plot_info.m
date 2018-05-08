%
% Author: Yosyp Schwab
%
% Function plots <x,y,z> data and optionally shows balls for atoms, sticks
% for bonds, and numbers of atomic indeces.
% Inputs:
%        bonds nearest neighbors connected graph <sq-matrix>
%        neighbs interatomic distances between all pairs <sq-matrix>
%        atom_graph with nodes = atom coords, bonds = edges <graph>
% Outputs:
%        fig1,fig2,fig3 figure handles
%
function [fig_bonds,fig_neighbs,fig_graph] = ...
    plot_info(bonds,neighbs, atom_graph)

    fig_bonds = figure;
        image(bonds,'CDataMapping','scaled');
        title('Adjacency Matrix of Bonds');
        xlabel('Atom #', 'FontWeight', 'bold', 'Color', 'black');
        ylabel('Atom #', 'FontWeight', 'bold', 'Color', 'black');
        xt = get(gca, 'XTick'); set(gca, 'FontSize', 16); set(gca, 'LineWidth', 2);
        
    fig_neighbs = figure;
        image(neighbs,'CDataMapping','scaled');   
        title('Nearest Neighbor Distance (darker=shorter)');
        xlabel('Atom #', 'FontWeight', 'bold', 'Color', 'black');
        ylabel('Atom #', 'FontWeight', 'bold', 'Color', 'black');
        xt = get(gca, 'XTick'); set(gca, 'FontSize', 16); set(gca, 'LineWidth', 2);        
    
    fig_graph = figure;
        plot(atom_graph);
        title('Undirected Graph of System');
        xt = get(gca, 'XTick'); set(gca, 'FontSize', 16); set(gca, 'LineWidth', 2);        
    highlight(atom_graph,24,'NodeColor','g','EdgeColor','g')

end
