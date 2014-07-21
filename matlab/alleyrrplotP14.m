% ALLEYRRPLOTP14 - plot Printer's Alley 140711 responses

% Created: 17-Jul-2014, JSA.
% Revised: 17-Jul-2014, JSA, v1.
% Version: v1.
% set(0, 'DefaultAxesFontSize', 16);
% set(0, 'DefaultAxesFontName', 'Arial');

%% initialization
% rpath = 'P140711/'; % response path, string
rpath = '~/Desktop/alley/audio/';
rname = 'alleyway140709_ir';   % response file name prefix, string
r = 6;  % response id, index

delta = 5;  % response plot duration, seconds

nbinst = 256;   % high time resolution analysis window length, nbins
nbinsf = 2048;  % high frequency resolution analysis window length, nbins
nskip = 64; % spectrogram hop size, samples

ptau = 80; % impulse response plot duration, milliseconds
peta = 68;  % impulse response start time, milliseconds


%% read impulse response
[ir1, fs] = wavread([rpath, rname, int2str(r), 'M1.wav']);
ir2 = wavread([rpath, rname, int2str(r), 'M2.wav']);

ntaps = round(delta*fs);
ir1 = ir1(1:ntaps);


%% plot impulse responses, spectrograms

% plot spectrograms
% figure(1);
% [~, ax1] = ftgram(ir1, fs, 'rir', 'nbins', nbinst, 'nskip', nskip, 'tanhbeta', 5, 'dbrange', 60);
% 
% figure(2);
% [~, ax2] = ftgram(ir1, fs, 'music', 'nbins', nbinsf, 'nskip', nskip);

% plot impulse response onsets
ptaps = round(ptau*fs/1000);
pstart = round(peta*fs/1000);
indexp = pstart + [0:ptaps-1];
figure(3);
% plot vertical lines where we expect reflections, according to our model
alleyirimage;
distances_of_interest = [3, 4, 7, 12, 1:4:56, 2:4:60]; 
% after reflection #10, we start to ignore reflections originating from the
% back of the speaker
h1 = plot((indexp-pstart)*1000/fs, 0.9*ir1(indexp)+1, 'b'); % grid;
hold all
h2 = plot((indexp-pstart)*1000/fs, ir2(indexp), 'color', [0.2 0.6 0.2]);
for i = distances_of_interest,
    plot(ones(1,51)*(distances(i)+0.4), -0.5:0.05:2.0, '-', 'color', [0.7 0.7 0.7]);
    hold all;
end;
% plot((indexp-pstart)*1000/fs, , 'color', [0.15 0.7 0.15]);
hold off;


xlabel('time, milliseconds'); ylabel('response | amplitude');
title('Impulse responses, Printers Ink Alley, orthogonal figure-8 mics');
% uistack(h2, 'top');
% uistack(h1, 'top');
legend('Mic 1', 'Mic 2');
