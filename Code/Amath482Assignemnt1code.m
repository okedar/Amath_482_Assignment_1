% Clean workspace 
clear all; close all; clc

%                   Data and Variables

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata
L = 15; % spatial domain 
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); %65 equaly spaced sections between -10 and 10
x = x2(1:n); 
y = x; 
z = x;
k = (2 * pi/(2 * L)) * [0:(n/2 - 1) -n/2:-1]; % assoasiates xy in 3d space to a frequency from 2pi
                                              % function to 2L function
ks = fftshift(k); 
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks,ks,ks); %frequency domain
ave = zeros(n,n,n); %added


%%
%                   Space-Time Signal Extraction

figure(1)
for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    M = max(abs(Un),[],'all');
    isosurface(X,Y,Z,abs(Un)/M,0.7) 
    axis([-20 20 -20 20 -20 20]) 
    grid on
    drawnow 
end

%%
%                   Averaging

figure(2)

for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Unt = fftn(Un);
    ave = ave + Unt; 
    M = max(abs(Unt),[],'all');
    isosurface(X,Y,Z,abs(Unt)/M,0.7) 
    axis([-20 20 -20 20 -20 20]) 
    grid on
    drawnow 
end

ave = ave/49;
ave = fftshift(abs(ave)); % shift matrix to accomidate shift 
aveMax = max(abs(ave),[],'all');

figure(3)
isosurface(Kx,Ky,Kz,abs(ave)/aveMax,0.7)
axis([-20 20 -20 20 -20 20]) 
grid on
drawnow

[Max_ave, Max_index] = max(abs(ave(:))); % Max_index gives index number of Max value in matrix ave
[max_kx, max_ky, max_kz] = ind2sub([n,n,n], Max_index);  % max_kx gives index location of on the kx
                                                         % at Max index (same me for ky and kz)

%%
%                   Filtering

tau = 0.2;
k0x = Kx(max_kx, max_ky, max_kz);  % Finds Kx value at Max index 
k0y = Ky(max_kx, max_ky, max_kz);  % Finds Ky value at Max index
k0z = Kz(max_kx, max_ky, max_kz);  % Finds Kz value at Max index
filter = fftshift(exp(-tau*(Kx - k0x).^2) .* exp(-tau*(Ky - k0y).^2) .* exp(-tau*(Kz - k0z).^2)); 

figure(4)
for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Unt = fftn(Un);
    Untf = filter.*Unt; 
    Unf = ifftn(Untf);
    [Max_ave, Max_index] = max(abs(Unf(:)));
    [max_x, max_y, max_z] = ind2sub([n,n,n], Max_index);
    location_x(j) = X(max_x, max_y, max_z);
    location_y(j) = Y(max_x, max_y, max_z);
    location_z(j) = Z(max_x, max_y, max_z);
end

figure(4)
plot3(location_x, location_y, location_z, 'Linewidth', 3);
title('Sub Trajectory');
xlabel('X Position');
ylabel('Y Position');
zlabel('Z Location'); 
XY_coordinates = [location_x.' , location_y.'];
grid on










