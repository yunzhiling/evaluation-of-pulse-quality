function [PPG,b,bb,rho1,rho2,rho3,c1,c2,c3]=evalpulse(data,time)
PPG=data;
fs=1/mode(diff(time));
% Quality test
[b,bb,rho1,rho2,rho3,PPG]=QPW(PPG,time,fs,[.5,18]);

ro1=[0,rho1(2:end-1),0];ro2=[0,rho2(2:end-1),0];ro3=[0,rho3(2:end-1),0];
% quality-based sorting
[a1,b1]=sort(ro1,'descend');
[a2,b2]=sort(ro2,'descend');
[a3,b3]=sort(ro3,'descend');
mbb=median(bb);mbb=round(fs*mbb);if mbb<=30;mbb=2*mbb;end
% Taking 5 best quality segments , using three different SQIs
for i=1:5
      [ma,mb]=findpeaks(PPG(floor(fs*(b(b1(i)))-10):ceil(fs*(b(b1(i)+1)))));
      if isempty(ma), [ma,mb]=findpeaks(PPG(floor(fs*(b(b1(i)))-10):ceil(fs*(b(b1(i)+2)))));end
      m(1,i)=mb(ma==max(ma))+floor(fs*(b(b1(i)))-10)-1;M1(:,i)=PPG(m(1,i)-20:m(1,i)+mbb-10);clear ma mb
      [ma,mb]=findpeaks(PPG(floor(fs*(b(b2(i)))-10):ceil(fs*(b(b2(i)+1)))));
       if isempty(ma), [ma,mb]=findpeaks(PPG(floor(fs*(b(b2(i)))-10):ceil(fs*(b(b2(i)+2)))));end
      m(2,i)=mb(ma==max(ma))+floor(fs*(b(b2(i)))-10)-1;M2(:,i)=PPG(m(2,i)-20:m(2,i)+mbb-10);clear ma mb
      [ma,mb]=findpeaks(PPG(floor(fs*(b(b3(i)))-10):ceil(fs*(b(b3(i)+1)))));
       if isempty(ma), [ma,mb]=findpeaks(PPG(floor(fs*(b(b3(i)))-10):ceil(fs*(b(b3(i)+2))))); end
      m(3,i)=mb(ma==max(ma))+floor(fs*(b(b3(i)))-10)-1;M3(:,i)=PPG(m(3,i)-20:m(3,i)+mbb-10);clear ma mb
%     [mm,m(1,i)]=max(PPG(floor(fs*(b(b1(i)))-10):ceil(fs*(b(b1(i)+1)))));m(1,i)=m(1,i)+floor(fs*(b(b1(i)))-10)-1;M1(:,i)=PPG(m(1,i)-20:m(1,i)+45);
%     [mm,m(2,i)]=max(PPG(floor(fs*(b(b2(i)))-10):ceil(fs*(b(b2(i)+1)))));m(2,i)=m(2,i)+floor(fs*(b(b2(i)))-10)-1;M2(:,i)=PPG(m(2,i)-20:m(2,i)+45);
%     [mm,m(3,i)]=max(PPG(floor(fs*(b(b3(i)))-10):ceil(fs*(b(b3(i)+1)))));m(3,i)=m(3,i)+floor(fs*(b(b3(i)))-10)-1;M3(:,i)=PPG(m(3,i)-20:m(3,i)+45);
end


% finding the averaged segment
c1=mean(M1')';cx1=min(find(diff(c1)>0));cX1=max(find(diff(c1)<0));
c2=mean(M2')';cx2=min(find(diff(c2)>0));cX2=max(find(diff(c2)<0));
c3=mean(M3')';cx3=min(find(diff(c3)>0));cX3=max(find(diff(c3)<0));


% % subplot(1,3,1),plot(c1(max([1,cx1-1]):cX1+1)),axis tight
% % subplot(1,3,2),plot(c2(max([1,cx2-1]):cX2+1)),axis tight
% % subplot(1,3,3),plot(c3(max([1,cx3-1]):cX3+1)),axis tight
%finding the features
if abs(fs-66.6667)>0.1
c1=resample(c1,0:1/fs:(length(c1)-1)/fs,1/0.015);
c2=resample(c2,0:1/fs:(length(c2)-1)/fs,1/0.015);
c3=resample(c3,0:1/fs:(length(c3)-1)/fs,1/0.015);
end
c1=c1(5:end-3);c2=c2(5:end-3);c3=c3(5:end-3);
figure

subplot(2,3,1),plot(diff(diff(c1)));hold on,plot([1,length(c1)],[0 0]),ylabel('APG'),title('averaged over 5 beats with highest SQI1'),axis tight
subplot(2,3,2),plot(diff(diff(c2)));hold on,plot([1,length(c2)],[0 0]),ylabel('APG'),title('averaged over 5 beats with highest SQI2'),axis tight
subplot(2,3,3),plot(diff(diff(c3)));hold on,plot([1,length(c3)],[0 0]),ylabel('APG'),title('averaged over 5 beats with highest SQI3'),axis tight
subplot(2,3,4),plot(c1),ylabel('PPG'),axis tight
subplot(2,3,5),plot(c2),ylabel('PPG'),axis tight
subplot(2,3,6),plot(c3),ylabel('PPG'),axis tight



