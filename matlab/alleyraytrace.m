% ALLEYRAYTRACE - simulate alleyway response via specular ray tracing.

% (c) Copyright 2014 Abel Innovations.  All rights reserved.
%
% Created: 18-May-2014, JSA - based on raytraceroom, v2.
% Revised: 16-July-2014, JSA, v1.
% Version: v1.


%% initialization

% system
pskip = 5; % samples skipped between frames, samples

c = 0.01/pskip;    % sound speed, meters per sample
nsamp = pskip*9600;    % total distance radiated, samples
nr = 4096;    % ray count, rays

% alley geometry
cant = 0.5;    % alley wall cant, meters
dip = 0.0;    % alley floor cant, meters
vertices = [cant 0 1 2 2-cant; 8 0 -dip 0 8];  % alley vertices, meters
nv = size(vertices,2);
alleytop = max(vertices(2,:));

% loudspeaker
speaker = [1/sqrt(5); exp(-0.5)];

% microphone
microphone = [sqrt(2); exp(-0.2)];
rho = 0.02; % microphone ray capture distance, meters

% film variables
framerate = 30; % film frame rate, frames per second

framesize = [720 540];  % film frame size, pixels
paperdpi = 150; % paper pixel density, dots per inch

rfscale = 1024; % rayfield resolution, pixels

% output variables
mflag = 1;  % plot microphone rays, indicator
writeflag = 1;  % write film frames, indicator
filmpath = 'raytrace_film'; % film frame directory, string
filmfile = [timecode, '.mp4'];  % output audio file name, string


%% form sides, side normals

% form sides
sides = [vertices(:,1:nv-1); vertices(:,2:nv)];

% compute normals, offsets
temp = sides(1:2,:) - sides(3:4,:);
normals = [temp(2,:); -temp(1,:)] ./ (ones(2,1) * sqrt(sum(temp.^2)));
offsets = sum(vertices(:,1:nv-1).*normals);


%% trace rays, plot results

% get figure
fid = figure(gcf);
set(fid, 'Position', [360 502 560*framesize/framesize(1)]);
set(fid, 'PaperUnits', 'inches', 'PaperOrientation', 'portrait');
set(fid, 'PaperPosition', [0 0 framesize/paperdpi]);
set(gca, 'Position', [0 0 1 1]);

% compute axis limits
temp = [min(vertices,[],2) max(vertices,[],2)];
axislimits = 0.2*kron([1 1], [-1 1]) + [temp(1,:) temp(2,:)];

% initialize rayfield parameters
vscale = max(vertices') - min(vertices');
voffset = min(vertices');
rfsize = ceil(rfscale*vscale/max(vscale));

hsize = floor(rfscale*rho/max(vscale));
h = hanning(2*hsize+1);
h = h*h'/sum(sum(h*h'));

% initialize ray
ray = speaker*ones(1,nr);

rayx = [ray(1,:); zeros(nsamp-1,nr)];
rayy = [ray(2,:); zeros(nsamp-1,nr)];

% initialize ray increment
theta0 = sort(2*pi*rand(1,nr));
delta = c*[cos(theta0); sin(theta0)];

deltax = [delta(1,:); zeros(nsamp-1,nr)];
deltay = [delta(2,:); zeros(nsamp-1,nr)];

% loop through samples
for t = [2:nsamp],
    
    % increment rays
    ray = ray + delta;

    % reflect rays, velocities
    for v = [1:nv-1],
        % test reflection
        projection = normals(:,v)'*ray - offsets(v);
        rflag = (projection < 0) .* (ray(2,:) < alleytop);

        % update ray
        ray = ray - 2*normals(:,v)*(projection .* rflag);

        % update delta
        delta = delta - 2*normals(:,v)*((normals(:,v)'*delta).*rflag);

    end;

    % assign rays, velocities
    rayx(t,:) = ray(1,:);
    rayy(t,:) = ray(2,:);

    deltax(t,:) = delta(1,:);
    deltay(t,:) = delta(2,:);

    % note captured rays
    lrflag = sum(abs(ray - microphone*ones(1,nr))) < rho;

    % plot geometry
    if ~rem(t,pskip),
        % set frame
        frame = t/pskip;

        % display progress
        if (rem(frame, 10*framerate) == 0),
            fprintf('%4d: %s\n', round(frame/framerate), time);
        elseif (rem(frame, framerate) == 0),
            fprintf('.');
        end;

        % plot frame
        if mflag,
            % include microphone paths
            plot(sides([1 3],:), sides([2 4],:), '.-r', ...
                speaker(1), speaker(2), '*g', ...
                ray(1,:), ray(2,:), '.b', rayx(1:t,lrflag), rayy(1:t,lrflag), '-c', ...
                microphone(1), microphone(2), 'og');

        else,
            % current ray points only
            plot(sides([1 3],:), sides([2 4],:), '.-r', ...
                speaker(1), speaker(2), '*g', ...
                ray(1,:), ray(2,:), '.b');

        end;
% %         axis(axislimits);
% %         axis('equal');
axis([-4.3252   6.3252  -0.2000  8.2000]);
        drawnow;

        % print frame
        if writeflag,
            print(gcf, '-djpeg', [filmpath, '/', num2str(frame, '%04d') '.jpg']);
        end;

    end;

end;


%% construct film

if writeflag,
    % generate film
    iv = [' -i ', filmpath, '/%04d.jpg'];
    r =  [' -r ', int2str(framerate)];
    system(['ffmpeg', r, iv, ' ', filmpath, '/', filmfile]);

    disp('');

end;

