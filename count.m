%
% Make CNT's
% http://turin.nss.udel.edu/research/tubegenonline.html
% post-process:
% cut header
% $ awk '{ print $2 " " $3 " " $4}' input.txt  > output.txt

clear;clc;
% CNT = dlmread('cnt.xyz', ' ');
% CNT = dlmread('longcnt.xyz', ' ');
CNT = dlmread('midcnt.xyz', ' ');
x = CNT(:,1);
y = CNT(:,2);
z = CNT(:,3);
natoms = length(x);

% marker size
s = 1*ones(1,length(x))';

% subplot(2,2,1);
hold on; grid on; view(40,35);
xlim([-5 4]);
ylim([-5 4]);
zlim([-5 4]);
scatter3(x,y,z,s,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]);
    
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
    nnindex = find(neighbs(i,:) == min(neighbs(i,:)));
    for k = nnindex
        line([x(i) x(k)], ...
             [y(i) y(k)], ...
             [z(i) z(k)], 'Color', 'red', 'LineWidth', 1);
        mytxt = text(x(i),y(i),z(i),sprintf('%d',i));         
        mytxt.Color = 'blue'; mytxt.FontSize = 16;
        bonds(i,k) = 1;
%          fprintf('Adding %d -> %d (%4.4f)\n', i , k, neighbs(i,k));         
    end
end
   
  
% subplot(2,2,2);
% image(bonds,'CDataMapping','scaled');
% 
% subplot(2,2,3);
% image(neighbs,'CDataMapping','scaled'); colorbar;

mygraph = graph(bonds);
% figure;plot(mygraph);

firstnode = 2;
% firstnode = 18;
% firstnode = 24;
N1 = firstnode;
N2 = NaN(3,2);
N3 = NaN(6,2);

N1 = neighbors(mygraph,firstnode);
k=1;
for i=1:length(N1)
    nb2 = neighbors(mygraph,N1(i));
    nb2 = nb2(nb2 ~= firstnode);
%     fprintf("\nN2: i:%d, length of nb2: %d\n", i, length(nb2));

    if length(nb2) < 2
        nb2 = [nb2; NaN];
    end    
    N2(i,:) = nb2;

    disp(N2(i,:))

    for j=1:length(N2(i,:))
        disp(N2(i,j));
        if isnan(N2(i,j)) == true
%             fprintf('true\n %d');
            N3(k,:) = [NaN NaN];
        else
            nb3 = neighbors(mygraph, N2(i,j));
            nb3 = nb3(nb3 ~= N1(i));  
%             fprintf("N3: k:%d, length of nb2: %d\n", k, length(nb3));   
%             disp(nb3)
            if length(nb3) < 2
                nb3 = [nb3; NaN];
            end
    %         disp(nb3)
            N3(k,:) = nb3;
        end
        k=k+1;
    end
    

end

rings = zeros(3,natoms);
bitrings = zeros(3,natoms);
count = 1;
% fprintf("Loop:c-i-j-k:\n");
for i=1:length(N1)
    for j=1:2
        for k=1:2
            if j == 1
                n3k = 2*i-1;
%                 disp(N3(n3k,k));
            else
                n3k = 2*i;
            end
%             if N3(n3k,k) > 0 && N2(i,j) > 0
%             if isnan(N3(n3k,k)) == false && isnan(N2(i,j)) == false
%             if isnan(N2(i,j)) == false && N2(i,j) > 0
%             if isnan(N3(n3k,k)) == false && N3(i,j) > 0
            if isnan(N3(n3k,k)) == false 
            if isnan(N2(i,j)) == false 
                fprintf("Loop:%d-%d-%d-%d: %d -> %d -> %d -> %d\n", ... 
                        count, i,j,k,firstnode, N1(i), N2(i,j), N3(n3k,k));
                bitrings(count,[firstnode, N3(n3k,k)]) = 1;
                rings(count,[firstnode, N1(i), N2(i,j) N3(n3k,k)]) = 1;
                count = count+1;    
            end
            end
        end
    end
end


twopair = nchoosek(1:count-1, 2);
thispair = [];
pairs = 0;
for i=1:length(twopair)
    fullrings(i) =  sum(bitand( ...
                            bitrings(twopair(i,1),:), ...
                            bitrings(twopair(i,2),:) ) );
    if fullrings(i) == 2
         thispair(end+1,:) = find(bitor( ...
                            rings(twopair(i,1),:), ...
                            rings(twopair(i,2),:) ));
        pairs = pairs+1;
    end
end

for i=1:pairs
    for j=1:6
        scatter3(x(thispair(i,j)),y(thispair(i,j)),z(thispair(i,j)),250,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1]);
    end
end
% image(rings,'CDataMapping','scaled');
% figure;image(fullrings,'CDataMapping','scaled');
% plot(fullrings);


disp(thispair);


