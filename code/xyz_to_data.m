clear;clc;

% CNT = dlmread('../input_files/cnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/longcnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/midcnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/cappedcnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/gencnt.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/c60.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/nanotube_carbon_lattice.xyz', ' ', 2, 1);
CNT = dlmread('../input_files/cnt-33-1110.xyz', ' ', 2, 1);
% CNT = dlmread('../input_files/1defect_cnt-33-1110.xyz', ' ', 2, 1);

% bond_length = 1.4210;
% bond_length = 1.334205;
bond_length = 1.2321/2;
x = CNT(:,1);
y = CNT(:,2);
z = CNT(:,3);
natoms = length(x);

x = x+abs(min(x));
y = y+abs(min(y));
z = z+abs(min(z));

% For carbon nanotube stretching only:
% Need to set the first few lattice widths at the top and bottom of the
% nanotube to rigid, so we need to identify those values here:
% rz(1): Rigid Z first (biggest)
% rz(2): Rigid Z second 
% rz(3): Rigid Z first (smallest)
% etc ...
rz(1) = max(z);
rz(2) = max(z(z<max(z)));
rz(3) = max(z(z<max(z(z<max(z)))));
rz(4) = min(z);
rz(5) = min(z(z>min(z)));
rz(6) = min(z(z>min(z(z>min(z)))));

fID = fopen('md-runs/test_cnt.data', 'w');

% Line 1: <#-of-atoms> <time>
fprintf(fID, '%d 0\n', natoms);

% Line 2: <system size x,y,z>
fprintf(fID, '%4.4f %4.4f %4.4f\n', ...
             max(x)+30*bond_length, ...
             max(y)+30*bond_length, ...
             max(z)+2*bond_length);

% Line 3: <system center>
fprintf(fID, '%4.4f %4.4f %4.4f\n', ...
             max(x)/2, ... 
             max(y)/2, ... 
             max(z)/2);

% Line 4: Type of each particle
for i=1:natoms
    fprintf(fID, '1 ');
end
fprintf(fID, '\n');

% Line 5: khist, x, y, z
for i=1:natoms
    if any(rz(:) == z(i))
        fprintf(fID, '3 %4.4f %4.4f %4.4f ', x(i), y(i), z(i));  
        fprintf('found rigid! %4.4f\n', z(i));
    else
        fprintf(fID, '0 %4.4f %4.4f %4.4f ', x(i), y(i), z(i));
    end
end
fprintf(fID, '\n');

% Line 6: velocity vx, vy, vz
for i=1:natoms
    fprintf(fID, '0 0 0 ');
end
fprintf(fID, '\n');

fclose(fID);
