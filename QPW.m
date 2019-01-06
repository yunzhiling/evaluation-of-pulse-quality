function [b,bb,rho1,rho2,rho3,pPPG]=QPW(PPG,t,fs,BK)
% Find the quality score of the pulsatile signal
% inputs:
% PPG, the pulsatile signal (vector)
% t, time (sec)
% fs, sampling frequency (Hz)
% range, an approximate range of valid signal (samples) for example [1:5000]
% outputs:
% b, detected beats (sec)
% bb, beat to beat intervals (sec)
% rho1, signal quality score by direct matching
% rho2, signal quality score by resampling
% rho3, signal quality score by DTW


% preprocessing

vt=find(t>0);
PPG=PPG(vt);
t=t(vt);
range=[fix(length(PPG)/4):fix(3*length(PPG)/4)];

% detecting the initial (transient) spikes..
ppg=abs(PPG);
as1=30*prctile(ppg,60)-29*median(ppg);
as2=30*prctile(ppg,30)-29*median(ppg);
PPG(ppg>as1)=mean(PPG(range));
PPG(ppg<as2)=mean(PPG(range));
%[fb,fa] = butter(2,[.5,7]/fs,'bandpass');
[fb,fa] = butter(2,BK/fs,'bandpass');
PPG = filtfilt(fb,fa,PPG);
%%% Finding the template (T) and beat onset (r)

[T,r]=Find_PW_Template(PPG,fs,range);
%%% Direct Matching:
for i=1:length(r)-2
   j=max(find(T.Beats(:,1)<=i));
   Tt=T.template{j};
   Rho1=corrcoef(Tt,PPG(r(i)+1:r(i)+length(Tt)));
   rho1(i)=Rho1(2,1);
   clear Tt
end
rho1(rho1<=0)=0;
%%% Linear sampling- stretch/compress
for i=1:length(r)-1
   j=max(find(T.Beats(:,1)<=i));
   Tt=T.template{j};
   if r(i+1)-r(i)<5*length(Tt);
   Rho2=corrcoef(Tt,resample(PPG(r(i)+1:r(i+1)),length(Tt),length(PPG(r(i)+1:r(i+1)))));
   rho2(i)=Rho2(2,1);
   else
       rho2(i)=0;       
   end
   clear Tt
end
rho2(rho2<=0)=0;
%% DTW
for i=1:length(r)-1
   j=max(find(T.Beats(:,1)<=i));
   Tt=T.template{j};
   if r(i+1)-r(i)<5*length(Tt);
   [Dist12,D12,k12,w]=dtw(Tt,PPG(r(i)+1:r(i+1))');
   dd=PPG(r(i)+1:r(i+1));
   DD(w(:,1))=dd(w(:,2));
   Rho3=corrcoef(Tt,DD);
   rho3(i)= Rho3(2,1);   
   clear DD dd Dist12 D12 k12 w
   else
       rho3(i)=0;       
   end
   clear Tt
end
rho3(rho3<=0)=0;
% %%% Clipping Detection
% for i=1:length(r)-1
%    rho41=sum(PPG(r(i)+1:r(i+1))==max(PPG(r(i)+1:r(i+1))))/(r(i+1)-r(i));
%    rho42=sum(PPG(r(i)+1:r(i+1))==min(PPG(r(i)+1:r(i+1))))/(r(i+1)-r(i));
%    rho4(i)=1-(rho41+rho42);
% end
% figure,hold on,
% plot(r(1:length(rho1)),100*rho1,'k'),
% plot(r(1:length(rho2)),100*rho2,'r'),
% plot(r(1:length(rho3)),100*rho3,'m'),
% plot(r(1:length(rho4)),100*rho4,'g'),
% figure,
% subplot(2,1,1),plot(PPG),hold on, plot(r,PPG(fix(r)),'r*'),
% subplot(2,1,2),plot(r(1:length(rho1)),100*rho1,'k'),hold on,
%plot(r(1:length(rho2)),100*rho2,'r'),
b=t(r(1:end-1));
bb=diff(b(b>0));
pPPG=PPG;    