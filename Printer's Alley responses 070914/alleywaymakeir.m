% alleymakeir - form alleyway impulse responses from 140709 measurements

% Created: 11-Jul-2014, JSA.
% Revised: 11-Jul-2014, JSA, v1.
% Version: 1.

%% initialization
sname = 'ssx_48_20.wav'; % sine sweep file name, string
[ss, fs] = wavread(sname);  % sine sweep, amplitude; sampling rate, Hz

rname = 'alleyway140709_20r';   % response file name prefix, string
nr = 13;    % measurement count, responses

ntaps = 10*fs;   % impulse response duration, samples

% measurement: speaker orientation, speaker position, mic pattern
%  1: N,  63N 20E, f8 f8
%  2: N,  63N 20E, f8 f8
%  3: E,  63N 20E, f8 f8
%  4: E,  63N 20E, f8 f8
%  5: W,  63N 20E, f8 f8
%  6: W,  63N 20E, f8 f8

%  7: N,  30N 20E, o f8
%  8: N,  30N 20E, o f8
%  9: N,  64N 20E, o f8
% 10: N,  64N 20E, o f8

% 11: N, 128N 20E, o f8
% 12: N,   0N 18W, o f8
% 13: N,   0N 18W, o f8


%% form impulse responses

% form sweep transform
sbins = 2^nextpow2(length(ss));
SS = repmat(fft(ss, sbins),1,2);

% initialize responses
ir1s = zeros(ntaps,nr);
ir2s = zeros(ntaps,nr);

% generate impulse responses
for r = [1:nr],
    % read response
    rr = wavread([rname, int2str(r), '.wav']);

    % form impulse response
    temp = real(ifft(fft(rr, sbins) ./ SS));
    ir1s(:,r) = temp(1:ntaps,1);
    ir2s(:,r) = temp(1:ntaps,2);

    % plot responses
    figure(r);
    ftgram([ir1s(:,r) ir2s(:,r)], fs, 'rir', 'tanhbeta', 5, 'nbins', 256, 'nskip', 32);

    drawnow;

end;


