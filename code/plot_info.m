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
        
    fig_neighbs = figure;
        image(neighbs,'CDataMapping','scaled');    
    
    fig_graph = figure;
        plot(atom_graph);

end
