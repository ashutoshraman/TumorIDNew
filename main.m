%Code to simulate Matt's laser moving across the surface of a sample
close all;
laser_power = 80; %0.200; %Amperes
spot_size = 3.0; %mm
myLaser = laser405(laser_power, spot_size); %myLaser is instance of laser405
myArgs.length = 40; %mm   % myArgs is a data struct array with length and width as args
myArgs.width = 40; %mm
mySample = sample405(myArgs); %mySample is instance of sample405 class, with myArgs struct
mySample.set_laser(myLaser); %mySample method of set_laser with myLaser instance as arg (class within class)

%set spectra
load('mean_healthy_spectra');
load('mean_tumor_spectra');
load('wavelengths');
mySample.tumor_spectra = mean_tumor_spectra;
mySample.healthy_spectra = mean_healthy_spectra;
mySample.wavelengths = wavelengths;
                                    
%Create hole for tumor and fill it w/ Tumor
length = max(size(mySample.x_coord));
width = max(size(mySample.y_coord));
thickness = mySample.thickness;

custom_lesion = create_lesion(mySample.state, length, width, thickness, 'circle');
mySample.tumor_state = custom_lesion;

% mySample.state(length/4:length/4 + length/2, width/4:width/4 + width/2) = thickness / 2; %carve out square of half of depth
% mySample.tumor_state = mySample.state;  % change this and above line to make circular lesion
I = find(mySample.tumor_state == thickness);
mySample.tumor_state(:,:) = thickness;
mySample.tumor_state(I) = NaN;

%Create a scan path
L = 0; %length (mm)
W = 30; %width (mm)
d = 1.75; %beam width (mm)
alpha = 1.0; %step size in terms of r between raster cuts
dir = 'x'; %direction of cuts 'x' or 'y';
speed = 1; %speed of cut mm/s
pp_mm = 1; %points per millimeter 
default_args = {L,W,d,alpha,dir,speed,pp_mm};
myScanPath = cut_path_obj('raster', default_args);
myScanPath.set_location([20,20]); %center of mySampleopen

% %Create a scan path
% myScanPath = cut_path_obj('empty', -1);
% myScanPath.set_x_points([5  10 20 20 28 32]);
% myScanPath.set_y_points([20 10 20 30 25 15]);
% myScanPath.p_points = [80 80 80 80 80 80];
% myScanPath.t_points = [0.1 0.2 0.3 0.4 0.5 0.6];
% %myScanPath.set_location([0,0]); %center of mySample

%Pre allocate space for acquired spectra
num_pts = max(size(myScanPath.t_points));
sz_wavelength = max(size(mySample.wavelengths));
mySample.acquired_spectra_series = repmat(struct('pts_xy',[NaN, NaN], 'spectra', nan(1,sz_wavelength)), num_pts, 1);                                    

%execute the path
mySample.perform_ablation(myScanPath); %not actually ablating, just same function which calls acquire spectra function

%Generate a laser spot
lspace = linspace(-1.5*myLaser.spot_size,1.5*myLaser.spot_size,99);
[X,Y] = meshgrid(lspace(1:2:end),lspace(1:2:end));
beam_profile = myLaser.pwr_function(sqrt((X.^2 + Y.^2)));
beam_profile(beam_profile < 0.01) = NaN;
beam_profile = beam_profile ./ max(max(beam_profile)); %normalize for plotting
%surf(beam_profile); axis square;




%Plot the results
H = figure(); %subplot(211);
plot_gcf(mySample, 'on', H); view(0,90);
hold on;
num_pts = max(size(myScanPath.x_points));
myScanPath.currentPoint = 1; %Reset count
k = 1;
while true
    next_pt = myScanPath.nextPoint; %[t, x, y]
    if isnan(next_pt), break; end
    %surf(X + next_pt(2), Y + next_pt(3), beam_profile + mySample.thickness);
    %shading interp
    scatter3(next_pt(2), next_pt(3), max(max(beam_profile)) + mySample.thickness, 2000, 'red', 'o');
    %scatter3(next_pt(2), next_pt(3), max(max(beam_profile)) + mySample.thickness, 10, 'red', 'o');
    text(next_pt(2) - .3, next_pt(3) + 5, max(max(beam_profile)) + mySample.thickness, num2str(k),'FontSize',15);
    k = k + 1;
end
axis equal
ylabel('Length (mm) [Y]','FontSize',15);
xlabel('Width (mm) [X]','FontSize',15);


%Smooth healthy and tumor spectra for comparison
smoothed_healthy_spectra = conv(mySample.healthy_spectra, [1 1 1 1 1 1 1 1] ./ 8);
smoothed_tumor_spectra = conv(mySample.tumor_spectra, [1 1 1 1 1 1 1 1] ./ 8);

% H = figure; %for the spectra
% for k = 1:num_pts
%     %subplot(ceil(num_pts/2), 2,k);
%     figure(H); H.clo;
%     plot(mySample.wavelengths, smoothed_healthy_spectra(8:end), 'k');
%     hold on;
%     plot(mySample.wavelengths, smoothed_tumor_spectra(8:end), 'g');
%     smoothed_spectra = conv(mySample.acquired_spectra_series(k).spectra, [1 1 1 1 1 1 1 1] ./ 8);
%     plot(mySample.wavelengths, smoothed_spectra(8:end), '--r');
%     title(sprintf('Spectra %d', k));
%     ylabel('Normalized Itensity A.U. (500nm -> 1.0)', 'FontSize', 15); 
%     legend({'Healthy','Tumor','Acquired'},'FontSize',15);    
%     xlabel('Wavelengths (nm)', 'FontSize', 15);
%     axis([0 1200 0 1.5]);
%     pause
% end

[Y, I_500] = min(abs((mySample.wavelengths - 500))); %I_x is index of wavelength nearest to x
[Y, I_545] = min(abs((mySample.wavelengths - 545))); %I_x is index of wavelength nearest to x
[Y, I_575] = min(abs((mySample.wavelengths - 575))); %I_x is index of wavelength nearest to x

%max ratio change between healthy and tumor between 475 and 700nm only
[Y, I_527] = min(abs((mySample.wavelengths - 527))); %I_x is index of wavelength nearest to x
[Y, I_687] = min(abs((mySample.wavelengths - 687))); %I_x is index of wavelength nearest to x

%max intensity change between healthy and tumor between 475 and 700nm only
[Y, I_527] = min(abs((mySample.wavelengths - 527))); %I_x is index of wavelength nearest to x
[Y, I_576] = min(abs((mySample.wavelengths - 576))); %I_x is index of wavelength nearest to x

%max difference 
[Y, I_471] = min(abs((mySample.wavelengths - 471))); %I_x is index of wavelength nearest to x
[Y, I_576] = min(abs((mySample.wavelengths - 576))); %I_x is index of wavelength nearest to x


mySpectra = mySample.acquired_spectra_series(1).spectra;

min1 = mySpectra(I_500) / mySpectra(I_545);
min2 = mySpectra(I_500) / mySpectra(I_575);
min3 = mySpectra(I_471) - mySpectra(I_576);

line1 = zeros(1,num_pts);
line2 = zeros(1,num_pts);
line3 = zeros(1,num_pts);
x_pos = zeros(1,num_pts);

for k = 1:num_pts
    %subplot(ceil(num_pts/2), 2,k);
    mySpectra = mySample.acquired_spectra_series(k).spectra;
    x_pos(k) = mySample.acquired_spectra_series(k).pts_xy(1);
    
    %scatter(x_pos, mySpectra(I_500) / mySpectra(I_545) - min1,'blue'); hold on; %make this more important
    line1(k) = mySpectra(I_500) / mySpectra(I_545) - min1;
    
    %scatter(x_pos, mySpectra(I_500) / mySpectra(I_575) - min2, 'red');
    line2(k) = mySpectra(I_500) / mySpectra(I_575) - min2;
    
    %scatter(x_pos, mean(mySpectra(I_425:I_750)),'green');
    %scatter(x_pos, mySpectra(I_527) / mySpectra(I_687), 'yellow');
    %scatter(x_pos, mySpectra(I_527)/mySpectra(I_687), 'magenta');
    
   
    %scatter(x_pos, mySpectra(I_471) - mySpectra(I_576) - min3, 'yellow');
    line3(k) = mySpectra(I_471) - mySpectra(I_576) - min3;
end

figure;
P1 = patch([0 10 10 0],[0 0 10 10], 'cyan'); set(P1,'facealpha',0.3); hold on;
P2 = patch([10 30 30 10],[0 0 10 10], 'magenta'); set(P2,'facealpha',0.3); hold on;

plot(x_pos, line1, 'blue'); hold on;
plot(x_pos, line2, 'green');
plot(x_pos, line3, 'red');

P3 = patch([30 40 40 30],[0 0 10 10], 'cyan'); set(P3,'facealpha',0.3); hold on;
%title('Tumor Edge Detection (center is tumor)','FontSize',15);
ylabel('Relative Intensity (a.u)', 'FontSize', 15); 
xlabel('X Position (mm)', 'FontSize', 15);
axis([0 40 0 0.5]);
legend({'Healthy','Tumor','500/545','500/575','471 - 576'},'FontSize',15)
axis square;

return
