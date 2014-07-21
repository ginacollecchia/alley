% ALLEYIRIMAGE - form alleyway impulse response using the image method.

% Created: 16-May-2014, JSA.
% Revised: 16-May-2014, JSA, v1.
% Version: v1.


%% initialization

% system
fs = 44100; % sampling rate, Hz
vp = 342.5;    % sound speed, meters per second (tuned to measured IR)

% alleyway geometry
cant = 0.04;    % alley wall cant, meters
dip = 0.0;    % alley floor cant, meters
vertices = [cant 0 2 2-cant; 8 0 0 8];  % alley vertices, meters
% vertices = [cant 0 1 2 2-cant; 8 0 -dip 0 8; 0 0 0 0 0];  % alley vertices, meters
nv = size(vertices,2);  % vertex count, vertices

gains = [1 1 1 0];    % surface reflectivities, fraction

% loudspeaker, microphone
% loudspeaker = [1.5936; 1.3462];  % source position, meters
% calling the left wall 0 on the x-axis, the ground 0 on the y-axis, and
% the source 0 on the z-axis, mic in front of it
loudspeaker = [1.3333; 1.2821]; % 0];
rho = 0.3;   % source radius, meters

% microphone = [1.5428; 0.7272]; % microphone position, meters
microphone = [1.1282; 1.6667; -1.6410];

% image method
order = 50; % maximum number of reflections computed per path, count
countmax = 2^20;    % maximum image count, images


%% generate images

% form sides, normals
sides = [vertices; [vertices(:,2:nv) vertices(:,1)]];

temp = sides(1:2,:) - sides(3:4,:);
normals = -[temp(2,:); -temp(1,:)] ./ (ones(2,1) * sqrt(sum(temp.^2)));
offsets = sum(vertices.*normals);

% form speaker images
temp = [1:nv];
indexk = temp(find(gains > 0));

images = [loudspeaker zeros(2,countmax-1)];
iside = zeros(1,countmax);
igain = [1 zeros(1,countmax)];
iparent = [1 zeros(1,countmax)];
visible = [1 zeros(1,countmax)];
count = 1;
s1 = 0;
for i = [1:order],
    % display progress
    if ~rem(i,5),
        fprintf('|');
    else
        fprintf('.');
    end;

    % loop through new images
	s0 = s1+1;
	s1 = count;
	for j = [s0:s1],
		image = images(:,j);

		% loop through polygon sides
		for k = indexk,
			normal = normals(:,k);
			rho = offsets(k);

			% form images
			d = rho - image'*normal;
			reflection = image + 2*d*normal;

			% validate image
			if (d >= 0),
				% valid image
				count = count + 1;
				images(:,count) = reflection;
				iparent(count) = j;
				iside(count) = k;
				igain(count) = gains(k) * igain(j);

				% check visibility
                visible(count) = 1;
				p = count;
				vmic = microphone(1:2);
				for n = [1:i],
					% find image, side
					vimage = images(:,p);
					vside = sides(:,iside(p));

					v1 = vside(1:2);
					v2 = vside(3:4);
                    pvp = eye(2) - (v2-v1)*(v2-v1)'/((v2-v1)'*(v2-v1));
                    c = vmic + (vimage-vmic)*(vimage-vmic)'*pvp*(v1-vmic)/((vimage-vmic)'*pvp*(vimage-vmic));

					% test visibility
					if ~((v1-v2)'*(v2-c) * (v1-v2)'*(v1-c) <= 0),
                    	visible(count) = 0;
					end;

					% reflect microphone
					normal = normals(:,iside(p));
					rho = offsets(iside(p));
					d = rho - vmic'*normal;
					vmic = vmic + 2*d*normal;

					% find reflection parent
					p = iparent(p);
				end;

			end;
		end;
	end;
end;
fprintf('\n');

% resize variables
igain = igain(1:count);
visible = visible(1:count);
images = images(:,1:count);

% form indices
indexi = find(~visible);
indexv = find(visible);

% plot results
figure(1);
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
% plot(vertices(1,[1:end 1]), vertices(2,[1:end 1]), '-k');
plot(vertices(1,:), vertices(2,:), '-k');
hold on;
plot(images(1,indexi), images(2,indexi), 'or', ...
    images(1,indexv), images(2,indexv), '*b'); grid;
plot(microphone(1), microphone(2), 'xk');
plot(loudspeaker(1), loudspeaker(2), '*k');
hold off;
title([num2str(cant,2),'m cant and ', num2str(dip,2),'m dip: ', int2str(sum(~visible)), ' invisible (red) and ', int2str(sum(visible)), ' visible (blue) reflections ']);
xlabel('distance, meters '); ylabel('distance, meters ');
axis([-40 40 -10 10]);




%% form impulse response

% compute, sort arrival distances
arrivals = [images(:,indexv) - microphone(1:2)*ones(1,length(indexv));  -microphone(3)*ones(1,length(indexv))];
distances = sqrt(sum(arrivals.^2))/vp * 1000; % print out in milliseconds
% images(:,indexv) % coordinates of visible images

[delays ndx] = sort(distances/vp);

% form arrival times, amplitudes
alpha = igain(ndx) ./ sqrt((1 + delays*vp/rho)/(1 + delays(1)*vp/rho));

