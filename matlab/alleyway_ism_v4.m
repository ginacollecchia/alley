clear all; close all;

%% setup
% all measurements in feet, so convert speed of sound to get samples later
% coordinate of source
x_s = 3.5/3.25; y_s = 5.5/3.25; h = y_s;
% geometry of alleyway: width and height
W = 6/3.25;
H = 24/3.25; 
% define maximum number of reflections to model
maxorder = 100;

%% define sources, count order to correctly apply filters
k = 1;
for order = -maxorder:maxorder-1,
    xcoord(k) = 2*order*W - x_s;
    xcoord(k+1) = xcoord(k) + 2*x_s;
    k = k + 2;
end;
% first two coordinates correlate to direct path and reflection off floor
xcoord = [xcoord, xcoord]; % xcoord repeated for reflections off floor
ycoord = y_s*[ones(1,maxorder*4), -ones(1,maxorder*4)];
% keep track of the number of reflections: this is the order for each row
% of virtual sources
orders = [maxorder:-1:1, maxorder:-1:1, 0, 0, 1:maxorder-1, 1:maxorder-1];

%% plot virtual sources
figure();
for i = 1:length(xcoord), % the length(xcoord)
    plot(xcoord, ycoord,'o');
    hold all;
end;

totL = (xcoord(3)-2):1:(xcoord(end)+2);
% plot the ground
plot(totL, zeros(1,length(totL)),'k', 'LineWidth',2);
% plot the alleyway walls
plot([0,0],[0,24/3.25],'k','LineWidth',2);
plot([W,W],[0,28/3.25],'k','LineWidth',2);

%% factor in a microphone
x_l = 3.5/3.25;
y_l = y_s;
% plot the source and mic on top of all of it so we can see em
plot(x_s,h,'ro','MarkerSize',5,'LineWidth',2);
plot(x_l,y_l,'rx','MarkerSize',5,'LineWidth',2);
grid on;
hold off;
title('Virtual Sources, Original Source, and Microphone in Alleyway');
% axis([xcoord(1)-2, xcoord(end)+2, ycoord(1)-2, ycoord(end)+2]);

% compute distance the sound travels from each source to the microphone
% coordinate
% dir_path_dist = sqrt((x_l-x_s)^2+(y_l-y_s)^2);
% dir_path_off_floor_dist = sqrt((x_l-x_s)^2+(y_l-(-y_s))^2);
% dist = [dir_path_dist, dir_path_off_floor_dist];

for i = 1:length(xcoord),
    dist(i) = sqrt((xcoord(i)-x_l)^2+(ycoord(i)-y_l)^2);
end;

% put direct path and its reflection off floor at the beginning of the
% matrix

%% sound absorption coeffs of painted concrete block (dense)
fs = 48000;
m = 20; % upsample to m*fs
nbins = 4096;
alpha = sqrt(1-[0.1, 0.05, 0.06, 0.07, 0.09, 0.08]); % calling 0.1 an "outlier"
fc = 125*2.^(0:5);
[b,a] = bilinear([alpha(5)/(2*pi*1000), alpha(2)],[1/(2*pi*1000),1],fs*m);
tfhat = freqz(b, a, pi*[0:nbins]'/nbins);
% plot the filter
% figure(); semilogx(fc,20*log10(alpha), 'o', [0:nbins]/nbins*fs*m/2, 20*log10(abs(tfhat)), '-'); grid;
% xlim([100 10000]);
% title('Absorption filter for painted dense concrete block');

%% reflection filter for ceiling
cv = 343; % meters per second
% fc_h = cv/W;
% [bh,ah] = butter(2,fc_h*2/(fs*m),'low');
% tfhat_h = freqz(bh, ah, pi*[0:nbins]'/nbins);
% figure(); semilogx(fc_h, 20*log10(0.5), 'o', [0:nbins]/nbins*fs*m/2, 20*log10(abs(tfhat_h)),'-'); grid;
% xlim([100, 10000]);

%% build the impulse response
implen = 2*round(max(max(dist))*fs*m/cv);
% store the impulse in imp
imp = zeros(implen,1);
diracsize = 48;
dirac = zeros(diracsize,1); dirac(1) = 1;
wallReflection = dirac * ones(1,max(orders)+1); % first column is impulse
% ceilingReflection = dirac * ones(1,max(orders)+1); % first column is impulse
k = 1; % diameter/aperture of balloon source

% precompute the filters; i is the order of the reflection
orders_h = orders;
for i = 1:max(orders), wallReflection(:,i+1) = filter(b,a,wallReflection(:,i)); end;

% for j = 1:max(orders_h), ceilingReflection(:,j+1) = filter(bh,ah,ceilingReflection(:,j)); end;

delay = round(dist*fs*m/cv);
% sfactor = 1./(1+dist/k).^(0.5);
sfactor = (1+dist/k).^(-0.1);

for nr = 1:length(orders),
    % for j = 1:length(orders),
        impWReflection = [zeros(1,delay(nr)), sfactor(nr) * wallReflection(:,orders(nr)+1)', ...
            zeros(1,implen-delay(nr)-length(dirac))]';
        % Add it to impulse response
        imp = imp + impWReflection;
    % end;
end;

% dist_h = ...

% and then add ceilingReflection
% for i = 1:length(orders_h),
%     for j = 1:length(orders_h),
%         impCReflection = [zeros(delay(i,j),1); sfactor(i,j) * ceilingReflection(:,orders_h(i)+1); zeros(implen-delay(i,j)-diracsize,1)];
%         % Add it to impulse response
%         imp = imp + impCReflection;
%     end;
% end;

%% convolve our impulse response imp with the recorded clap
[clap,fsclap] = wavread('clap.wav');
clap = resample(clap, fsclap, fs);
% downsample imp after normalizing
imp = imp./(0.9*max(abs(imp)));
imp_48k = resample(imp, 1, m);
y = fftfilt([zeros(480,1);clap],imp_48k); 
y = y./(0.9*max(abs(y)));
% y = y(1:fs);
figure(); 
ftgram(y,fs,'rir');
title('Synthesized impulse response');

