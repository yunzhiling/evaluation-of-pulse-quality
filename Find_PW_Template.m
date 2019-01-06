function [Te,r]=Find_PW_Template(PPG,fs,range)
% Find the template of the pulsatile signal
% inputs:
% PPG, the pulsatile signal (vector)
% fs, sampling frequency (Hz)
% range, an approximate range of valid signal (samples) for example [1:5000]
% outputs:
% Te, templates
% r, detected beats (sample)

% transpose if necessary
if size(PPG,1)>=size(PPG,2),PPG=PPG';end
%100,1
%
%normalize
Kk=range; 
Kk=Kk(abs(PPG(Kk))>-1);
PPG=(PPG-mean(PPG(Kk)))./std(PPG(Kk));
%
% detect beat onset
r = wabp(PPG', 0,.005,fs);r=round((fs/125)*r);
%figure,plot(PPG),hold on, plot(r,PPG(fix(r)),'r*'),title('PPG'),
%
%figure
Temp1_invalid=0;
first_valid_Temp=1;
%%% Initial template generator
for i=1:ceil(length(PPG)/1500)
%     AR(i,:)=autocorr(PPG(1500*(i-1)+1:min(1500*i,length(PPG)))',min(length(PPG(1500*(i-1)+1:min(1500*i,length(PPG))))/2,250));% find the autocorrelation of each 30 sec window
%     [AR_p{i},AR_peaks{i}]=findpeaks(AR(i,:));
%     AR_Diff{i}=diff(AR_peaks{i});
%    L(i)=AR_Diff{i}(1);% find the average beat length in each 30 sec window
% find the range of beats
rm=1500*(i-1)+1;
rM=min(1500*i,length(PPG));
if length(r(r>=rm & r<rM))>1.5
L(i)=round(mean(diff(r(r>=rm & r<rM))));
if i==ceil(length(PPG)/1500),L(i)=round(min(diff(r(r>=rm & r<rM))));end
    r_min=sum(r<=(1500*(i-1)))+1;
    %r_max=sum(r<=(1500*i));
    r_max=sum(r<=(1500*i-L(i)+1));
    Te.Beats(i,:)=[r_min,r_max];
    tm=zeros(r_max-r_min+1,L(i));
    for j=1:r_max-r_min
        tm(j,:)=PPG(r(r_min+j-1):r(r_min+j-1)+L(i)-1);
    end 
    if size(tm,1)>1.5,T1{i}=mean(tm);else T1{i}=tm;end % The first template (T1): average of all the beats in the 30s window with each beat beginning at the onset of the beat and ending at the length of the template.
    %figure,title('Initial templates'),plot(T1{i}),hold on
    %%% Correlation of template 1 with each beat
    for j=1:r_max-r_min
        C=corrcoef(T1{i}',PPG(r(r_min+j-1):r(r_min+j-1)+L(i)-1));
        c(i,j)=C(1,2);
        clear C
    end
    if size(tm,1)<1.5
        Valid_ratio(i)=0;
    else 
    if size(c,1)==i
        if sum(c(i,:)>0)>2
    %%% removing uncorrelated beats and recalculating the template
    % more than 0.6 of beats removed??
            Valid_ratio(i)=(sum(c(i,:)>=0.8))/(r_max-r_min);
        else 
            Valid_ratio(i)=0;
        end
    else
        Valid_ratio(i)=0;
    end
    end
    else
        Valid_ratio(i)=0;
    end
    if Valid_ratio(i)>=0.4
    T{i}=mean(tm((c(i,:)>=0.8)',:));
    end
    if Valid_ratio(i)<0.4 && i<1.5 
        Temp1_invalid=1;
    end
    if Temp1_invalid==0 && Valid_ratio(i)<0.4
        T{i}=T{i-1};
    end
    if Temp1_invalid==1 && i>1.5 && Valid_ratio(i)>0.4
        T{i-1}=T{i};
        Temp1_invalid=0;
        first_valid_Temp=i;
    end
clear tm  
end
for i=1:first_valid_Temp-2
    T{i}=T{first_valid_Temp};
end
%figure,for i=1:length(T),plot(T{i}),hold on,end
Te.template=T;
