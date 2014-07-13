clear all; close all;

%% setup
% all measurements in feet, so convert speed of sound to get samples later
x_s = 1.5;
y_s = 3;
h = y_s;
W = 6;
H = 24; 

maxorder = 100; % well acutally this is half the maxorder
orders = [];
% compute the coordinates of the virtual sources: first x-coordinate
xcoord = []; ycoord = [];
theta = pi*0.554/180; % canting angle of the left wall

%%  define sources, count order to correctly apply filters
k = 1;
for i = 1:maxorder, % go up to 6*4 sources, plus the first one at -x_s (so, 25)
    % wall = 2*i*W;
    % first to the left, first to the right, second on the left, second on
    % the right
%     sources_x = [-2*(i-1)*W-x_s, wall-x_s, -wall+x_s, wall+x_s];
%     xcoord = [xcoord, sources_x];
    % CANTING
    wall = 2*i*W*cos(theta);
    % first to the left, first to the right, second on the left, second on
    % the right
    sources_x = [-(wall-2*W*cos(theta))-x_s, wall-x_s, -wall+x_s, wall+x_s];
    xcoord = [xcoord, sources_x];
    sources_y = [y_s+k*sin(theta), y_s-k*sin(theta), y_s+(k+1)*sin(theta), y_s-(k+1)*sin(theta)];
    ycoord = [ycoord, sources_y];
    orders = [orders, k, k, k+1, k+1];
    k = k+2;
end;

% duplicate the matrix to get the xcoordinates of sources reflecting off
% ground (same xcoordinates)
xcoord = [xcoord, xcoord]; % since i left off the first order reflection
% we know all the heights: height of the source and -height of the source
% ycoord = y_s*[ones(1,maxorder*4), -ones(1,maxorder*4)];
% CANTING:
ycoord = [ycoord, -ycoord];
% even though yes, these do go through the floor (another order), we 
% assume perfect reflection off floor
orders = [orders, orders];

%% now factor in a microphone
x_l = 3.5;
y_l = y_s;

% include direct path and first reflection off of the floor
xcoord = [x_s, x_s, xcoord];
ycoord = [y_s, -y_s, ycoord];
orders = [0, 0, orders]; % again assuming perfect reflection off floor

% compute distance the sound travels from each source to the microphone
% coordinate
for i = 1:length(xcoord),
    dist(i) = sqrt((xcoord(i)-x_l)^2+(ycoord(i)-y_l)^2);
    % plot([xcoord(i),x_l],[ycoord(i),y_l],'k--');
end;

dist = dist/3.25; % convert to meters

%% plot virtual sources + mic/source
figure();
plot(xcoord,ycoord,'bo','MarkerSize',5);
hold all;
% plot the source reflecting off floor
plot(x_s,-h,'go','MarkerSize',5); % not added in yet
totL = (sources_x(3)-2):1:(sources_x(4)+2);
% plot the ground
plot(totL, zeros(1,length(totL)),'k', 'LineWidth',2);
% plot the alleyway walls
plot([0,0],[0,24],'k','LineWidth',2);
plot([W,W],[0,30],'k','LineWidth',2);
axis([sources_x(3)-2, sources_x(4)+2, -h-2, h+2]);
% plot the source and mic on top of all of it so we can see em
plot(x_s,h,'ro','MarkerSize',5,'LineWidth',2); % still not factored in yet to dist
plot(x_l,y_l,'rx','MarkerSize',5,'LineWidth',2);
grid on;
xlabel('Feet'); ylabel('Feet');
hold off;
title('Virtual Sources, Original Source, and Microphone in Alleyway');

%% now build the impulse response
fs = 48000;
m = 20;
fs_big = m*fs;
nbins = 4096;
cv = 343; % meters per second

%% sound absorption coeffs of painted concrete block (dense)
alpha = sqrt(1-[0.1, 0.05, 0.06, 0.07, 0.09, 0.08]);
fc = 125*2.^(0:5);
% [b,a] = bilinear([alpha(end)/1000/2/pi/1.05, alpha(1)],[1/1000/2/pi,1],fs);
[b,a] = bilinear([alpha(5)/(2*pi*1000), alpha(2)],[1/(2*pi*1000),1],fs*m);
tfhat = freqz(b, a, pi*[0:nbins]'/nbins);
figure(); semilogx(fc,20*log10(alpha), 'o', [0:nbins]/nbins*fs*m/2, 20*log10(abs(tfhat)), '-'); grid;
xlim([100 10000]);
title('Absorption filter for painted dense concrete block');

%% build the impulse response
implen = round(max(dist)/cv*fs);

imp = zeros(m*implen,1);
dirac = zeros(48,1); dirac(1) = 1;
wallReflection = dirac * ones(1,max(orders)+1); % first column is impulse
k = 1; % diameter/aperture of balloon source

% precompute the filters
for i = 1:max(orders), wallReflection(:,i+1) = filter(b,a,wallReflection(:,i)); end;

delay = round(dist*fs*m/cv);
sfactor = (1+dist/k).^(-0.1);

for nr = 1:length(orders),
    impReflection = [zeros(1,delay(nr)), sfactor(nr) * wallReflection(:,orders(nr)+1)', ...
        zeros(1,length(imp)-delay(nr)-length(wallReflection(:,orders(nr)+1)))]';
    % Add it to impulse response
    impReflection = impReflection(1:length(imp));
    imp = imp + impReflection;
end;

%% now add "ceiling" reflection
fc_h = cv/W; % cutoff frequency is related to the width of the alleyway...
% determine how many vertical reflections will happen in our model given
% the size
maxorder_h = round(max(dist)*3.25/(2*H)); % convert dist back to feet
k = 1;
for i = 1:maxorder_h,
    ycoord_h(k) = 2*i*H-y_s-y_l; % coordinates from above
    ycoord_h(k+1) = -2*i*H-y_s-y_l; % coordinates from below
    orders_h(k) = i;
    orders_h(k+1) = orders_h(k);
    k = k+2;
end;
xcoord_h = x_s*ones(1, length(ycoord_h));
% now compute distance for "ceiling" reflections
for i = 1:length(xcoord_h),
    dist_h(i) = sqrt((xcoord_h(i)-x_l)^2 + (ycoord_h(i)-y_l)^2);
end;

dist_h = dist_h/3.25; % convert to meters
implen_h = round(max(dist_h)/cv * fs);
imp_h = zeros(m*implen,1);
dirac_h = zeros(48,1); dirac_h(1) = 1;
% ceilingReflection = dirac_h * ones(1,max(orders_h)+1); % first column is impulse
k = 1; % diameter/aperture of balloon source

% [bh,ah] = butter(2,fc_h*2/fs);
% precompute the filters for ceiling reflections
% for i = 1:max(orders_h), ceilingReflection(:,i+1) = filter(bh,ah,ceilingReflection(:,i)); end;
ceilingReflection = dirac_h;

delay_h = round(dist_h*fs*m/cv);
sfactor_h = (1+dist_h/k).^(-0.9); % 0.9 is hand-wavy, just as the 0.1 is before

for nr = 1:length(orders_h),
    impReflection_h = [zeros(1,delay_h(nr)), sfactor_h(nr) * ceilingReflection',... % ceilingReflection(:,orders_h(nr)+1)', ...
        zeros(1,length(imp_h)-delay_h(nr)-length(ceilingReflection))]';
    % Add it to impulse response
    impReflection_h = impReflection_h(1:length(imp_h));
    imp_h = imp_h + impReflection_h;
end;

imp = imp + imp_h;

%% convolve with clap from measurement
[clap,fsclap] = wavread('clap.wav');
clap = resample(clap, fsclap, fs);
y = resample(imp, 1, m);
yy = fftfilt([zeros(480,1);clap],y); 
yy = yy./max(abs(yy));
% yy = yy(1:fs);
figure();
% plot([0:length(yy)-1]/(fs), yy);
ftgram(yy, fs, 'rir');
% grid on;
xlabel('Time (s)');
title('Synthesized impulse response');

% %% convolve with clap from measurement
% y_h = resample(imp_h, 1, m);
% yy_h = fftfilt([zeros(480,1);clap],y_h); 
% yy_h = yy_h./max(abs(yy_h));
% % yy_h = yy(1:fs);
% figure();
% % plot([0:length(yy)-1]/(fs), yy);
% ftgram(yy_h, fs, 'rir');
% % grid on;
% xlabel('Time (s)');
% title('Synthesized impulse response: ceiling only');
