% Alleyway IR analysis
% clear all;

%% geometry for ISM
% of alleyway
Lx = 30.0; % 30 mteters long
Ly = 4.0; % 4 meters wide
Lz = 10.0; % 10 meters high
numAlleys = length(Lx);

micPosX = 0.0; % mic in the middle both length-wise and width-wise
micPosY = 0.0; 
micPosZ = 2.0-Lz/2; % 2 meters off the ground
numMic = length(micPosX);

sourcePosX = 0.2;
sourcePosY = 0.0;
sourcePosZ = 2.0-Lz/2;
dirSor = [sin(180*pi/180) cos(180*pi/180) 0];

% get our T60s
% [x,Fs] = wavread('~/Documents/MATLAB/alley/balloon2_rightchan.wav');
% xr = resample(x,48000,Fs);
% fc = [200,400,800,1600,3200,6400,12800]; % 6 different bands
% fT60 = zeros(1,length(fc)-1);
% for i = 1:length(fc)-1,
%     fT60(i) = sqrt(fc(i)*fc(i+1)); % center frequencies
% end;
% [m,mindex] = max(abs(xr));
% % figure(1);
% T60 = calcT60(xr,48000,fc,mindex,100,12)'; % T60 is 6x1
fT60 = [0.001,400,800,1600,3200,6400,16000,50000]';
T60 = [15, 8, 4, 1.8, 1.4, 1.2, 0.5, 0.16]';
numT60 = size(T60,2);

%% filtering and sample rates
fclean = [100 7200];

%% 3-D model
% lengthFilt = 2.^(nextpow2(round(max(T60)*fs))-1); % in orig code it is only 1/3 second
lengthFilt = 8*8192;
fs = 16000;
c = 343;
distLimit = 4*c; % limits the impulse response (max(T60) = 11.5, too big)
data = zeros(numMic, numT60, numAlleys, lengthFilt);

for alley = 1:numAlleys,
    N0 = ceil(distLimit./[Lx(alley) Ly(alley) Lz(alley)]);
    N = N0(1); M = N0(2); K = N0(3);
    % max number to use in the simulation. Note N, M and K need to be different integers
    NTot = (2*N+1)*(2*M+1)*(2*K+1); % total number of rooms - used to size vectors

    V = (Lx(alley)*Ly(alley)*Lz(alley)); % volume of room
    
    % dmax=sqrt((Lx*(N+1))^2+(Ly*(M+1))^2+(Lz*(K+1))^2);%calculating the maximum observer to source distance (and hence max delay time)
    % tmax=dmax/c; % maximum time delay
    
    % distance between mics and source. Used to ensure filter lengths are long enough
    sor2micDist = sqrt((micPosX(alley)-sourcePosX(alley)).^2+(micPosY(alley)-sourcePosY(alley)).^2+(micPosX(alley)-sourcePosZ(alley)).^2); 
    numSamp = ceil(fs*distLimit/c)+ceil(fs*sor2micDist/c); % number of samples that will contain all the reflections

    micRotateAng = 0; %angle in radians
    rotMat = [cos(micRotateAng) sin(micRotateAng); 
        -sin(micRotateAng) cos(micRotateAng)];

    %% defining the list of image rooms
    count = 0;
    RM = zeros(NTot,3);
    for n = -N:N,
        for m = -M:M,
            for k = -K:K,
                count = count+1;
                RM(count,:) = [n,m,k];
            end;
        end;
    end;

    % RM = RM(find(sqrt(RM.^2*[Lx(room) Ly(room) Lz(room)].^2')<distLimit),:);
    % figure(1)
    % plot3(Lx(room)*RM(:,1),Ly(room)*RM(:,2),Lz(room)*RM(:,3),'+')

    %% information for calculating the sub sample filter
    len = 16;
    FSS = [0:len-1]*fs/len;
    [B, A] = butter(2,0.9);
    LP = freqz(B,A,FSS,fs);

    % loop to calculate the reflection filters defining the frequencies and T60s
    numTap = 15;
    ns = 1024;
    F0 = [0:ns-1]'*fs/ns;
    Delta = 1/fs;

    % for generic T60 investigations
    for t60 = 1:numT60,
        wT60 = 2*pi*fT60;
        % reflection coefficients
        HD = 10.^(-3./(343*T60(:,t60)./(4*V/(Lx(alley)*Ly(alley)*2+Lx(alley)*Lz(alley)*2+Lz(alley)*Ly(alley)*2)))); 

        % interpolating the reflection coefficients for all frequencies
        Des = interp1(log(wT60),HD,log(2*pi*F0(1:ns/2)),'linear');
        Des(1) = Des(2);

        T = [ones(ns/2,1), 2*cos(2*pi*F0(1:ns/2)*[1:numTap]*Delta)];
        hb = pinv(T)*Des;
        hbD = [flipud(hb(2:numTap+1)); hb];
        Hopt = freqz(hbD,1,F0,fs);
        figure(2342);
        semilogx(F0,(abs(Hopt)),F0(1:ns/2),(abs(Des)));
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        grid on;
        axis([0 20000 0 1]);

        ref = hbD';
        % ref = exp(-[0:9]*0.25)./sum(exp(-[0:9]*0.25))*0.95;
        Ref = ref; % defines the reflection filter
        maxLenFilt = 1+(length(ref)-1)*(N+M+K+1);
        filt = zeros(K+N+M+1,maxLenFilt);
        filt(1,1) = 1;
        for i = 2:K+N+M+1,
            filt(i,1:length(Ref)) = Ref;
            Ref = conv(Ref,ref);
        end;

        %% looping through source and mic locations
        for mic = 1:numMic
            posOb = [[micPosX(mic), micPosY(mic)]*rotMat, micPosZ(mic)]; % location of mics
            posSor = [sourcePosX(mic), sourcePosY(mic), sourcePosZ(mic)]; % location of source


            % calculating the locations of all of the images
            figure(2);
            POSSor = [RM(:,1)*Lx(alley)+(-1).^RM(:,1)*posSor(1), RM(:,2)*Ly(alley)+(-1).^RM(:,2)*posSor(2), RM(:,3)*Lz(alley)+(-1).^RM(:,3)*posSor(3)];%locations of all image sources
            hold all;
            plot3(POSSor(:,1),POSSor(:,2),POSSor(:,3),'o');
            plot3(posOb(1),posOb(2),posOb(3),'*');
            hold off;

            % calculating the distance to the observer
            sor2Ob = [posOb(1)-POSSor(:,1), posOb(2)-POSSor(:,2), posOb(3)-POSSor(:,3)]; % vector from sources to observer
            d = sqrt(sum(sor2Ob.^2,2)); % distance to the observer

            sampDelay = (d/c*fs); % equivalent time delay in samples

            % making impulse response vector
            imp = zeros(numSamp+maxLenFilt,1);

            % narrow down a set number of impulses
            nearSet = find((d-min(d)) < distLimit);

            % loop to go through all of the reflections (can't think of a way to do it
            % without a loop)
            for counter = 1:length(nearSet),
                i = nearSet(counter);
                DirSor = [(-1)^RM(i,1)*dirSor(1), (-1)^RM(i,2)*dirSor(2), (-1)^RM(i,3)*dirSor(3)]; % orientation of the image source relative to the original room
                DirectivityFactor = 0.5*(sum(sor2Ob(i,:)./sqrt(sum(sor2Ob(i,:).^2)).*DirSor)+1); % making a cardioid
                % DirectivityFactor=1;

                % dealing with sub sample delay
                % determines where to add the impulse in time
                st = round(sampDelay(i)); 
                % residual sub sample delay
                delay = st-sampDelay(i); 
                ph = exp(1i*delay*2*pi*FSS./fs).*LP;
                ph2 = [ph(1:len/2), 0, conj(fliplr(ph(2:len/2)))];
                SSDelay = circshift(ifft(ph2.'),5);

                % defining impulse
                imageImp = filter(SSDelay,1,filt(abs(RM(i,1))+abs(RM(i,2))+abs(RM(i,3))+1,:)');
                % putting together the directivity, the decay and the absorption
                imp(st:st+maxLenFilt-1) = imp(st:st+maxLenFilt-1)+(DirectivityFactor/d(i))*imageImp; 
            end;
            
            figure(213);
            plot(imp);
            title(['Mic = ' num2str(mic) ',  t60 = ' num2str(t60) ',  alley = ' num2str(alley)]);
            drawnow
           % pause

            data(mic,t60,alley,1:length(imp))=imp;
        end
    end
end

%% filtering the impulse responses
[BBP,ABP] = butter(2,fclean*2/fs);%bandpass filter to apply to the microphone filters
han = hanning(800);
win = [ones(lengthFilt-400,1);han(401:800)];
si = size(data,4);
if si > lengthFilt,
    data = data(:,:,:,1:lengthFilt);
    si = lengthFilt;
end;
for mic = 1:numMic,
    for t60 = 1:numT60, 
        for alley = 1:numAlleys, 
            temp = filter(BBP,ABP,circshift(squeeze(data(mic,t60,alley,:)), 30-find(data(mic,t60,alley,:) < 0,1)));
            IR(:,mic,t60,alley) = [temp ;zeros(lengthFilt-si,1)] .*win;
        end
    end
end

figure();
[z,~] = ftgram(IR,16000,'rir');
