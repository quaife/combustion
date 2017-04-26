clear
clc
close all

% loading files

set_data =1;

switch set_data
    case 1
        load('snap_4451.mat')
        L =10:10:870;
    case 2
        load('snap_2451.mat')
        L =10:10:490;
    case 3
        load('snap_1030.mat')
        L =10:10:320;
end


data_number=size(snap,2);
dt=.1;
t=0:dt:dt*data_number;

count = 0;

% DMD approximation

for j = L
    [dmdbasis y0 omega Atilde A] = dmd_comp_Q(snap(:,1:end-1),snap(:,2:end),j,dt);
    
    for i=1:data_number
        sol_dmd(:,i)=y0.*exp(omega*t(i));
    end
    count = count + 1
    %keyboard
    sol_dmd_full = real(dmdbasis*sol_dmd);
    error(count) = norm(sol_dmd_full-snap)/norm(snap);
    archivo = strcat('Reduced_solution_',num2str(set_data),'_',num2str(j),'.txt');
    save(archivo,'sol_dmd_full','-ascii')
    clear sol_dmd sol_dmd_full
end

% save DATA

%archivo = strcat('Reduced_solution_',num2str(set_data),'.txt');
%save(archivo,'sol_dmd_full')
archivo = strcat('Error_',num2str(set_data),'.txt');
save(archivo,'error')