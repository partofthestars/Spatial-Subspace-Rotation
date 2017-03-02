%%########################################################################
%%
%%	SpatialSubspaceRotation
%%
%%
%% Author   : Christian S. Pilz
%% Date     : 24.02.2017
%%
%% Contact  : cpi@partofthestars.com
%% Web Page : www.partofthestars.com
%%
%% Version  : Alpha RA 1.0
%%
%% License  : GPL v3
%%
%%########################################################################
%%
%%	subspaceRotation.m:
%%
%% Description:
%%
%% A Novel Algorithm for Remote Photoplethysmography: Spatial Subspace Rotation
%%
%% W. Wang, S. Stuijk and G. de Haan, 
%% "A Novel Algorithm for Remote Photoplethysmography: Spatial Subspace Rotation,"
%% in IEEE Transactions on Biomedical Engineering, vol. 63, no. 9, pp. 1974-1984, Sept. 2016.
%%
%%

clear all; close all;

base_system_dir='./data/images/';

load('./data/CMS50E_ppg.mat');
load('./data/diffusion_process.mat');

% analysis loop
%

frames=1452;

for f=1:frames
    
    file=[base_system_dir num2str(f+24) '.png'];
    
    I = imread(file);
    [rows cols dim]=size(I);
    
    [out bin]=generate_skinmap(file);
    
    r=double(I(:,:,1));
    r=r(bin==1);
    
    g=double(I(:,:,2));
    g=g(bin==1);
    
    b=double(I(:,:,3));
    b=b(bin==1);

    values=[r g b];
    
    %raw mean traces
    trace(f,:)=mean(values);
    
    %spatial RGB correlation:
    C=(values'*values)/(rows*cols);
    [V,D] = eig(C);
    
    [U_,S_,V_] = svd(C);
    
    V_tmp=V;

    [D,I] = sort(diag(D),'descend');
     V = V(:, I);

    U{f}=V;
    Sigmas{f}=D;
 
    if f>1
        %rotation between the skin vector and orthonormal plane
        R{f-1}=[U{f}(:,1)'*U{f-1}(:,2) U{f}(:,1)'*U{f-1}(:,3)];
         
        %scale change
        S{f-1}=[sqrt(Sigmas{f}(1)/Sigmas{f-1}(2)) sqrt(Sigmas{f}(1)/Sigmas{f-1}(3)) ];
        SR{f-1}=S{f-1}.*R{f-1};
        SR_backprojected{f-1}=SR{f-1}*[U{f-1}(:,2)'; U{f-1}(:,3)'];
        
        signal_u(f-1,:)=U{f}(:,1)';
        signal_r(f-1,:)=R{f-1};
        signal_s(f-1,:)=S{f-1};
        signal_sr(f-1,:)=SR{f-1};
        signal_sr_b(f-1,:)=SR_backprojected{f-1};
    end
    
end

 sigma=std(signal_sr_b(:,1))/std(signal_sr_b(:,2));  
 p=signal_sr_b(:,2);%-sigma*signal_sr_b(:,3);
 pp=signal_sr_b(:,1)-sigma*signal_sr_b(:,2);
 pp=pp-mean(pp);
 
 fs=25;
 windowSize=3;%sec
 overLap=2;%sec
 
 [blocks_one]=buffer(signal_sr_b(:,1),windowSize*fs,windowSize*fs-1)';
 [blocks_two]=buffer(signal_sr_b(:,2),windowSize*fs,windowSize*fs-1)';
 [blocks_three]=buffer(signal_sr_b(:,3),windowSize*fs,windowSize*fs-1)';
 
 [frames dim]=size(blocks_one);

for i=1:frames
    sigma=std(blocks_one(i,:))/std(blocks_two(i,:));  
	
    p_block=blocks_one(i,:)-sigma*blocks_two(i,:);
    pulse(i,:)=double(p_block-mean(p_block));
end


%low-pass filter
%
cutoff=2.0;%hz
fNorm = cutoff / (fs/2);                    
[b,a] = butter(10, fNorm, 'low');        
                                        
p_blocked = filtfilt(b, a, double(pulse(:,end)));

raw=trace(:,2);
raw = filtfilt(b, a, double(raw(:,end)));
raw = diff(raw);

%pulse oximeter
%

fs=60;

cutoff=2.0;%hz
fNorm = cutoff / (fs/2);                    
[b,a] = butter(10, fNorm, 'low');
ppg = filtfilt(b, a, double(ppg(:,end)));
cutoff=0.5;%hz
fNorm = cutoff / (fs/2);
[b,a] = butter(10, fNorm, 'high');                                      
ppg = filtfilt(b, a, double(ppg(:,end)));
ppg = diff(ppg);

%plot results
%

figure;

fs=25;

y=raw';
L=length(y);
NFFT = 2^nextpow2(L);
f = fs/2*linspace(0,1,NFFT/2+1);

subplot(4,2,1);
plot(raw);
title('Derivatives of green-channel means');
xlabel('Frames');
subplot(4,2,2);
spectrogram(raw,128,120,f,fs,'yaxis');
title('Spectrogram')
ylim([0 3]);

subplot(4,2,3);
plot(p_blocked);
title('Spatial Subspace Rotation')
xlabel('Frames');

y=p_blocked';
L=length(y);
NFFT = 2^nextpow2(L);
f = fs/2*linspace(0,1,NFFT/2+1);

subplot(4,2,4);
spectrogram(p_blocked,128,120,f,fs,'yaxis');
title('Spectrogram')
ylim([0 3]);

subplot(4,2,5);
plot(p_blocked_dp);
title('Diffusion Process')
xlabel('Frames');

y=p_blocked_dp;
L=length(y);
NFFT = 2^nextpow2(L);
f = fs/2*linspace(0,1,NFFT/2+1);

subplot(4,2,6);
spectrogram(p_blocked_dp,128,120,f,fs,'yaxis');
title('Spectrogram')
ylim([0 3]);

subplot(4,2,7);
plot(ppg);
title('CMS50E Pulse Oximeter')
xlabel('Frames');

fs=60;
y=ppg;
L=length(y);
NFFT = 2^nextpow2(L);
f = fs/2*linspace(0,1,NFFT/2+1);

subplot(4,2,8);
spectrogram(ppg,128,120,f,fs,'yaxis');
title('Spectrogram')
ylim([0 3]);
