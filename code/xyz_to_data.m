clear;clc;

% CNT = dlmread('../input_files/cnt.xyz', ' ');
% CNT = dlmread('../input_files/longcnt.xyz', ' ');
% CNT = dlmread('../input_files/midcnt.xyz', ' ');
% CNT = dlmread('../input_files/cappedcnt.xyz', ' ');
% CNT = dlmread('../input_files/gencnt.xyz', ' ');
CNT = dlmread('../input_files/c60.xyz', ' ');
x = CNT(:,1);
y = CNT(:,2);
z = CNT(:,3);
natoms = length(x);

fID = fopen('cnt.data', 'w');

% Line 1: <#-of-atoms> <time>
fprintf(fID, '%d 0\n', natoms);

% Line 2: <system size x,y,z>
fprintf(fID, '%4.4f %4.4f %4.4f\n', ...
        min(x) + max(x), ...
        min(y) + max(y), ...
        min(z) + max(z));

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
    fprintf(fID, '0 ');
end
fprintf(fID, '\n');

fclose(fID);
