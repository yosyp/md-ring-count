clear;clc;

% CNT = dlmread('../input_files/cnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/longcnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/midcnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cappedcnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/gencnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/c60.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/nanotube_carbon_lattice.xyz', ' ', 2, 1);
CNT = dlmread('../input_files/cnt-33-1110.xyz', ' ', 2, 1);

x = CNT(:,1);
y = CNT(:,2);
z = CNT(:,3);
natoms = length(x);

fID = fopen('md-runs/cnt.data', 'w');

% Line 1: <#-of-atoms> <time>
fprintf(fID, '%d 0\n', natoms);

% Line 2: <system size x,y,z>
fprintf(fID, '%4.4f %4.4f %4.4f\n', ...
        abs(min(x)) + max(x), ...
        abs(min(y)) + max(y), ...
        abs(min(z)) + max(z));

% Line 3: <system center>
fprintf(fID, '0 0 0\n');

% Line 4: Type of each particle
for i=1:natoms
    fprintf(fID, '1 ');
end
fprintf(fID, '\n');

% Line 5: khist, x, y, z
for i=1:natoms
    fprintf(fID, '0 %4.4f %4.4f %4.4f ', x(i), y(i), z(i));
end
fprintf(fID, '\n');

% Line 6: velocity vx, vy, vz
for i=1:natoms
    fprintf(fID, '0 0 0 ');
end
fprintf(fID, '\n');

fclose(fID);
