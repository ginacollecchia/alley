clear all; % close all;

%% setup
% all measurements in feet, so convert speed of sound to get samples later
x_s = 1.5;
y_s = 4;
h = y_s;
L = 6;
% Reflections
% reflections = struct('R',S,'order',0,'walls',[]);

maxorder = 10; % well acutally this is half the maxorder
orders = [];
% compute the coordinates of the virtual sources: first x-coordinate
xcoord = []; 
% we know all the heights: height of the source and -height of the source
ycoord = y_s*[ones(1,maxorder*4), -ones(1,maxorder*4)];

%%  define sources, count order to correctly apply filters
k = 1;
for i = 1:maxorder, % go up to 6*4 sources, plus the first one at -x_s (so, 25)
    wall = 2*i*L;
    sources = [(wall-2*L)-x_s, wall-x_s, -wall+x_s, wall+x_s];
    xcoord = [xcoord, sources];
    orders = [orders, k, k, k+1, k+1];
    k = k+2;
end;

% duplicate the matrix to get the xcoordinates of sources reflecting off
% ground (same xcoordinates)
xcoord = [xcoord, xcoord];

% even though yes, these do go through the floor (another order), we 
% assume perfect reflection off floor
orders = [orders, orders];

%% plot virtual sources
figure();
plot(xcoord,ycoord,'bo','MarkerSize',5);
hold all;
% plot the source reflecting off floor
plot(x_s,-h,'go','MarkerSize',5); % not added in yet
totL = (sources(3)-2):1:(sources(4)+2);
% plot the ground
plot(totL, zeros(1,length(totL)),'k', 'LineWidth',2);
% plot the alleyway walls
plot([0,0],[0,h+2],'k','LineWidth',2);
plot([L,L],[0,h+2],'k','LineWidth',2);
axis([sources(3)-2, sources(4)+2, -h-2, h+2]);


%% now factor in a microphone
x_l = 3.5;
y_l = 4.2;

% compute distance the sound travels from each source to the microphone
% coordinate
for i = 1:length(xcoord),
    dist(i) = sqrt((xcoord(i)-x_l)^2+(ycoord(i)-y_l)^2);
    plot([xcoord(i),x_l],[ycoord(i),y_l],'k--');
end;

% plot the source and mic on top of all of it so we can see em
plot(x_s,h,'ro','MarkerSize',5,'LineWidth',2); % still not factored in yet to dist
plot(x_l,y_l,'rx','MarkerSize',5,'LineWidth',2);

grid on;
hold off;
title('Virtual Sources, Original Source, and Microphone in Alleyway');

%% now build the impulse response
fs = 48000;
cv = 343; % meters per second

%% sound absorption coeffs of painted concrete block (dense)
alpha = sqrt(1-[0.1, 0.05, 0.06, 0.07, 0.09, 0.08]);
fc = 125*2.^(0:5);
[b,a] = bilinear([alpha(end)/1000/2/pi/1.05, alpha(1)],[1/1000/2/pi,1],fs);
figure();
freqz(b,a,4096,fs)
hold all
plot(fc,20*log10(alpha),'r+');
hold off;
title('Absorption filter for painted dense concrete block');

% compute direct path
dir_path_dist = sqrt((x_l-x_s)^2+(y_l-y_s)^2);
dir_path_off_floor_dist = sqrt((x_l-x_s)^2+(y_l-(-y_s))^2);
dist = [dir_path_dist, dir_path_off_floor_dist, dist];
dist = dist/3.25; % convert to meters
orders = [0, 0, orders]; % again assuming perfect reflection off floor

%% build the impulse response
implen = 4*round(max(dist)/cv*fs);

imp = zeros(implen,1);
dirac = imp; dirac(1) = 1;
for nr = 1:length(orders),
    delay = round(dist(nr)/cv*fs);
    impReflection = dirac;
    % apply the reflection filters the order of the reflection many times
    if orders(nr)>0,
        for i = 1:orders(nr), 
            impReflection = filter(b,a,impReflection);
        end;
    end;
    % Apply propogation delay: does natural e^-delay enveloping
    impReflection = filter([zeros(1,delay),1],1,impReflection);
    % Apply energy decay: this is decaying far too fast
    % impReflection = impReflection/dist(nr);
    % Add it to impulse response
    imp = imp + (-1)^(nr-1)*impReflection;
end;
figure();
plot([0:length(imp)-1]/fs, imp);
grid on;
xlabel('Time (s)');
title('Synthesized impulse response');
