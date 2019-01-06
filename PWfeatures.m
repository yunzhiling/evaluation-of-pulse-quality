function [F,T]=PWfeatures(c,mbb)
%%subplot(3,1,1),plot(c),hold on,xlim([0,length(c)]),title('PW')
[p1,p2]=findpeaks(c);
%%subplot(3,1,1),plot(p2,p1,'r*');
[v1,v2]=findpeaks(-c);
%%subplot(3,1,1),plot(v2,-v1,'ko');
for i=2:length(c)-1
    D1(i)=0.5*(c(i+1)-c(i-1));
end
D1(1)=D1(2);D1(end+1)=D1(end);
[pd1,pd2]=findpeaks(D1);
[vd1,vd2]=findpeaks(-D1);
zD1=[];
for i=2:length(D1)-1
    if D1(i-1)>0 & D1(i+1)<0
        zD1=[zD1,i];
    end
end
fzd1=find(diff(zD1)<1.5);
zD1(fzd1+1)=-1;
zD1=zD1(zD1>0.5);
for i=2:length(D1)-1
    D2(i)=0.5*(D1(i+1)-D1(i-1));
end
D2(1)=D2(2);D2(end+1)=D2(end);
[pd21,pd22]=findpeaks(D2);
[vd21,vd22]=findpeaks(-D2);
%%subplot(3,1,2),plot(D1),hold on,xlim([0,length(D1)]),title('1st derivative of PW')
%%subplot(3,1,2),plot(pd2,pd1,'o'),
%%subplot(3,1,2),plot(vd2,-vd1,'r*'),
%%subplot(3,1,2),plot(zD1,D1(zD1),'ko'),
%%subplot(3,1,3),plot(D2),hold on,xlim([0,length(D2)]),title('APW, 2nd derivative of PW')
%%subplot(3,1,3),plot(pd22,pd21,'o'),
%%subplot(3,1,3),plot(vd22,-vd21,'r*'),
f={p1,p2,-v1,v2,pd1,pd2,-vd1,vd2,zD1,pd21,pd22,-vd21,vd22};
% systolic peak amplitude - preceding valley amplitude
sysP=[max(p1),p2(max(find(p1>=(max(p1)))))];sysV=[-v1(max(find(v2<sysP(2)))),v2(max(find(v2<sysP(2))))];
EndPPG=max(find(c<=min(c(sysP(2):end))));
if abs(sysP(2)-EndPPG)<7
sysP=[p1(1),p2(1)];sysV=[-v1(max(find(v2<sysP(2)))),v2(max(find(v2<sysP(2))))];    
end
if abs(sysP(2)-EndPPG)<10
EndPPG=length(c);
end
if isempty(sysV)
    sysV=[c(1),1];
end
F(1)=sysP(1)-sysV(1);
F(2)=sysP(2)-sysV(2);
% Pulse width (half peak amplitude)
cx=c-sysV(1);
N=find(cx<=0.5*(sysP(1)-sysV(1)));
ph1=max(N(find(N<sysP(2))));
ph1=max(ph1,2);
x1=ph1-1;
x2=ph1+1;
y1=cx(ph1-1);
y2=cx(ph1+1);
a=0.5*(sysP(1)-sysV(1));
PH1=x1+(a-y1)*(x2-x1)/(y2-y1);
ph2=min([min(N(find(N>sysP(2)))),length(cx)-1]);
x1=ph2-1;x2=ph2+1;y1=cx(ph2-1);y2=cx(ph2+1);a=0.5*(sysP(1)-sysV(1));
PH2=x1+(a-y1)*(x2-x1)/(y2-y1);
F(3)=PH2-PH1;

% area


F(4)=trapz(c(sysV(2):EndPPG)-min(c(sysV(2):EndPPG)));

% notch
[nnh,ntch]=findpeaks(D2(sysP(2):EndPPG));ntch=ntch+sysP(2)-1;
Me=mean([sysP(2),EndPPG])-4;
Dsn=abs(ntch-Me);
    [Kn1,Kn2]=min(Dsn);
    Ntch=ntch(Kn2);NT=Ntch;
if ~isempty(v2(v2>sysP(2)&v2<EndPPG))
    DsN=abs(v2((v2>sysP(2)&v2<EndPPG))-Me);
    KN=v2((v2>sysP(2)&v2<EndPPG));
    [KN1,KN2]=min(DsN);
    Ntch=KN(KN2);
end
if isempty(Ntch)
    if length(v2)>1.5
    Ntch=v2(2);
    else
        Ntch=NaN;
    end
end


% 2nd drvtve

a1=NaN;a2=NaN;
b1=NaN;b2=NaN;
c1=NaN;c2=NaN;
d1=NaN;d2=NaN;
e1=NaN;e2=NaN;

[a1,aa2]=max(pd21(pd22<sysP(2))); a2=pd22(aa2);if isempty(a2),a2=pd22(1);end
b2=min(vd22(find(vd22>a2)));b1=D2(b2);
e1=D2(NT);e2=NT;
if isempty(a1),a1=NaN;a2=NaN;end
if isempty(e1),e1=NaN;e2=NaN;end
if isempty(b1),b1=NaN;b2=NaN;end
if sum(pd22>b2 & pd22<(e2-1))>0.5
    cc2=pd22(pd22>b2 & pd22<(e2-1));
    c2=cc2(1);c1=D2(c2);
    if sum(vd22>c2 & vd22<(e2-1))>0.5
        dd2=vd22(vd22>c2 & vd22<(e2-1));
        d2=dd2(1);d1=D2(d2);
    end
else
    if sum(vd22>b2 & vd22<(e2-1))>0.5
        dd2=vd22(vd22>b2 & vd22<(e2-1));
        d2=dd2(1);d1=D2(d2);
    end  
end

if ~(abs(Ntch)>-1)
    Ntch=e2;
end

% area ratio
if (abs(Ntch)>-1)
F(5)=trapz(c(Ntch:EndPPG)-min(c(sysV(2):EndPPG)))/trapz(c(sysV(2):Ntch-1)-min(c(sysV(2):EndPPG)));
else
    F(5)=NaN;
end
% beat-beat and pulse intervals and ratios
F(6)=mbb;
F(7)=-sysV(2)+EndPPG;
F(8)=F(2)/F(7);
% Augmontation index
if (abs(Ntch)>-1)
dp=vd22(vd22>sum(Ntch)+0 & vd22<EndPPG);
if ~isempty(dp)
    DP=dp(1)-1;
else
    [ddp,DP]=min(D2(Ntch:end));DP=DP+Ntch-2;
end
dp2=p2(p2>sum(Ntch)+0 & p2<EndPPG);
if ~isempty(dp2)
    DP=dp2(1);
end
%%plot(c),hold on,plot(DP,c(DP),'ro'),plot(Ntch,c(Ntch),'kx'),
F(9)=c(DP)-sysV(1);
else
    F(9)=NaN;
end
F(10)=F(9)/F(1);
F(11)=(F(1)-F(9))/F(1);
if (abs(Ntch)>-1)
F(12)=DP-sysV(2);
F(13)=Ntch-sysV(2);
F(14)=c(Ntch)-sysV(1);
F(15)=DP-sysP(2);
else
F(12)=NaN;
F(13)=NaN;
F(14)=NaN;
F(15)=NaN; 
end


% 
% plot(c),hold on,plot(D2,'r'),
% % plot(sysP(2),sysP(1),'x')
% % plot(sysV(2),sysV(1),'x')
% plot(a2,a1,'x')
% plot(b2,b1,'o')
% plot(c2,c1,'kx')
% plot(d2,d1,'ko')
% plot(e2,e1,'r*')

F(16)=b1/a1;
F(17)=c1/a1;
F(18)=d1/a1;
F(19)=e1/a1;
F(20)=(b1-c1-d1-e1)/a1;
F(21)=(b1-e1)/a1;
F(22)=(b1-c1-d1)/a1;
F(23)=(c1+d1-b1)/a1;
F(24)=b2-a2;
F(25)=e2-b2;

t1=abs(b1)/3;t2=abs(b1)/4;t3=abs(b1)/8;t4=abs(b1)/4;t5=abs(b1)/8;

if (a1>0) && (b1<0) && (d1<0) && (e1>0)
    if c1>0
        T=1;
    elseif d1-b1>t1 
        T=2;
    elseif d1-b1<=t1 && d1-b1>=t3 && c1-b1>=t2 
        T=3;
    elseif d1-b1<=t1 && d1-b1>=t3 && c1-b1<t2 
        T=4;
    elseif d1-b1<t3 && c1-b1>=t4 
        T=5;
    elseif d1-b1<t3 && c1-b1<t4 && c1-b1>=t5
        T=6;
    elseif d1-b1<t3 && c1-b1<t5
        T=7;
    elseif ~(abs(c1)>-1)
        T=7;        
    end    
else
    T=0;
end

