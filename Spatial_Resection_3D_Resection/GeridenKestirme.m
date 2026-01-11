clc ; clear all ; close all; format long g ; delete *.asv;
% Uzay geriden kestirme ile bir foto'nun dış yöneltme parametreleri
% bulunur.(omega,fi,kappa,XL,YL,ZL)
% Gerekenler : XYZ nesne uzayı koordinatları bilinen en az 3 kontrol noktası
%% Verilenler

%Distorsiyon
f=152.916;
disp(['f= ',num2str(f),' mm'])
f=f.*10^-3;

%Okuma
photo = xlsread('Koordinatlar_Uzay_Geriden_Kestirme.xlsx','B:C'); 
Xf=photo(:,1).*10^-3;Yf=photo(:,2).*10^-3;

XYZ = xlsread('Koordinatlar_Uzay_Geriden_Kestirme.xlsx','D:F');  
X=XYZ(:,1);Y=XYZ(:,2); Z=XYZ(:,3);
ng=size(XYZ,1);  % matris boyutu

%Deltalar bBaşlangıs
% dengelemeye sokacağım delta değerleri
delta1=[1 1 1 1];
delta2=[1 1 1 1];

%Omega Fi Başlangıç
% Dikkat : Soruda DÜŞEYE YAKIN BİR HAVA FOTOĞRAFI denildiği için..
% omega ve fi'yi aşağıda 0 aldık.
om(1)=0;
fi(1)=0;


%% ZL Hesabı
%Ho=[649.5450;640.6980;648.7437;650.2681;657.2923;648.7251];
a=1;  %delta 1 BI Yavlışl

for i=1:1:ng;
    for j=1:1:ng;
        if i==j || j<i;
            continue
        else
            AB2=(X(j,1)-X(i,1))^2+(Y(j,1)-Y(i,1))^2;
            syms H;
            Formul=(( (Xf(j,1)/f)*(H-Z(j,1)) - (Xf(i,1)/f)*(H-Z(i,1)) )^2 + ( (Yf(j,1)/f)*(H-Z(j,1)) - (Yf(i,1)/f)*(H-Z(i,1)) )^2)-AB2;
            Ho(a)=double(solve(Formul,H>0));
            a=a+1;
        end
    end
end
ZL=mean(Ho);
% ZL = izdüşüm merkezinin koordinatı
% Ho ' da tum yukseklikler var. Burada 6 tane var.
% Cunku soruda 4 nokta var. 4'ün 2'lisi = 6
disp(['Ortalama Yükseklik (ZL): ', num2str(ZL)]);

%% AX ve AY
for i=1:1:ng;
    AX(i,1)=(Xf(i,1)/f)*(ZL-Z(i,1)); 
    AY(i,1)=(Yf(i,1)/f)*(ZL-Z(i,1)); 
end

%% A matrisi
a=1;
for i=1:1:ng;

    A(a,1)=AX(i,1);
    A(a+1,1)=AY(i,1);

    A(a,2)=-AY(i,1); 
    A(a+1,2)=AX(i,1);

    A(a,3)=1;
    A(a+1,3)=0;

    A(a,4)=0;
    A(a+1,4)=1;
    
    a=a+2;
end

 a=1;   
 %% x ve v hesabı ve kapa
 for i=1:1:ng;

    l(a,1)=X(i,1); 
    l(a+1,1)=Y(i,1);
    a=a+2;

 end
 
x=inv(A'*A)*(A'*l); 
v=A*x-l;  

a=x(1,1);
b=x(2,1); 
XL=x(3,1); 
YL=x(4,1);

kap(1)=wrapTo2Pi(atan2(b,a))*180/pi

%% İterasyon
ii=0;

while max(abs(delta1))>0.0000005 || max(abs(delta2))>0.0001
    ii=ii+1;
    
    disp([num2str(ii),' nci iterasyon...'])
 
 

% Donusluk Matrisi Dönsturme Aşamaları:

 m(1,1)=cosd(fi(ii))*cosd(kap(ii));
 m(1,2)=sind(om(ii))*sind(fi(ii))*cosd(kap(ii))+cosd(om(ii))*sind(kap(ii));
 m(1,3)=-cosd(om(ii))*sind(fi(ii))*cosd(kap(ii))+sind(om(ii))*sind(kap(ii));
 m(2,1)=-cosd(fi(ii))*sind(kap(ii));
 m(2,2)=-sind(om(ii))*sind(fi(ii))*sind(kap(ii))+cosd(om(ii))*cosd(kap(ii));
 m(2,3)=cosd(om(ii))*sind(fi(ii))*sind(kap(ii))+sind(om(ii))*sind(kap(ii));
 m(3,1)=sind(fi(ii));
 m(3,2)=-sind(om(ii))*cosd(fi(ii));
 m(3,3)=cosd(om(ii))*cosd(fi(ii));
 
% Değer atamaları
 
 a=1;
 
 for i=1:1:ng;
 
    r(i,1)=m(1,1)*(X(i,1)-XL(ii))+m(1,2)*(Y(i,1)-YL(ii))+m(1,3)*(Z(i,1)-ZL(ii));
    s(i,1)=m(2,1)*(X(i,1)-XL(ii))+m(2,2)*(Y(i,1)-YL(ii))+m(2,3)*(Z(i,1)-ZL(ii));% 10 ninde 2
    q(i,1)=m(3,1)*(X(i,1)-XL(ii))+m(3,2)*(Y(i,1)-YL(ii))+m(3,3)*(Z(i,1)-ZL(ii));% 10 ninde 2
    
    b(a,1)=f/q(i,1)^2*(r(i,1)*(-m(3,3)*(Y(i,1)-YL(ii))+m(3,2)*(Z(i,1)-ZL(ii)))-q(i,1)*(-m(1,3)*(Y(i,1)-YL(ii))+m(1,2)*(Z(i,1)-ZL(ii))));
    
    b(a,2)=f/q(i,1)^2*(r(i,1)*(cosd(fi(ii))*(X(i,1)-XL(ii)) + sind(om(ii))*sind(fi(ii))*(Y(i,1)-YL(ii)) - cosd(om(ii))*sind(fi(ii))*(Z(i,1)-ZL(ii))) - q(i,1)*(-sind(fi(ii))*cosd(kap(ii))*(X(i,1)-XL(ii)) + sind(om(ii))*cosd(fi(ii))*cosd(kap(ii))*(Y(i,1)-YL(ii)) - cosd(om(ii))*cosd(fi(ii))*cosd(kap(ii))*(Z(i,1)-ZL(ii))));
    
    b(a,3)=-f/q(i,1)*(m(2,1)*(X(i,1)-XL(ii))+m(2,2)*(Y(i,1)-YL(ii))+m(2,3)*(Z(i,1)));

    b(a,4)=f/q(i,1)^2*(r(i,1)*m(3,1)-q(i,1)*m(1,1));

    b(a,5)=f/q(i,1)^2*(r(i,1)*m(3,2)-q(i,1)*m(1,2));

    b(a,6)=f/q(i,1)^2*(r(i,1)*m(3,3)-q(i,1)*m(1,3));
    
    b(a+1,1) = f/q(i,1)^2*(s(i,1)*(-m(3,3)*(Y(i,1)-YL(ii))+m(3,2)*(Z(i,1)-ZL(ii)))-q(i,1)*(-m(2,3)*(Y(i,1)-YL(ii))+m(2,2)*(Z(i,1)-ZL(ii)))); 
    
    b(a+1,2)=f/q(i,1)^2*(s(i,1)*(cosd(fi(ii))*(X(i,1)-XL(ii))+sind(om(ii))*sind(fi(ii))*(Y(i,1)-YL(ii))-cosd(om(ii))*sind(fi(ii))*(Z(i,1)-ZL(ii)))-q(i,1)*(sind(fi(ii))*sind(kap(ii))*(X(i,1)-XL(ii))-sind(om(ii))*cosd(fi(ii))*sind(kap(ii))*(Y(i,1)-YL(ii))+cosd(om(ii))*cosd(fi(ii))*sind(kap(ii))*(Z(i,1)-ZL(ii))));
   
    b(a+1,3)=f/q(i,1)*(m(1,1)*(X(i,1)-XL(ii))+m(1,2)*(Y(i,1)-YL(ii))+m(1,3)*(Z(i,1)-ZL(ii))); 
    b(a+1,4)=f/q(i,1)^2*(s(i,1)*m(3,1)-q(i,1)*m(2,1));
    
    b(a+1,5)=f/q(i,1)^2*(s(i,1)*m(3,2)-q(i,1)*m(2,2));

    b(a+1,6)=f/q(i,1)^2*(s(i,1)*m(3,3)-q(i,1)*m(2,3));
    

    J(i,1)=Xf(i,1)+f*r(i,1)/q(i,1); %Değerde hata
    K(i,1)=Yf(i,1)+f*s(i,1)/q(i,1);
    
    E(a,1)=J(i,1);
    E(a+1,1)=K(i,1);

    a=a+2;


 end

%% Yeni Omega Fi Kapa
B=[b(:,1) b(:,2) b(:,3) -1.*b(:,4) -1.*b(:,5) -1.*b(:,6) ];

Oxx=inv(B'*B);

BI=Oxx*(B'*E); 
V=B*BI-E;

delta1=BI(1:3);
delta2=BI(4:6);

D(:,ii)=BI(1:3)*(180/pi)*3600;

F(:,ii)=BI(4:6);

%Sonuçlar  İterasyon sebebiyle atamalar


om(ii+1)=om(ii)+BI(1)*180/pi;

fi(ii+1)=fi(ii)+BI(2)*180/pi;

kap(ii+1)=kap(ii)+BI(3)*180/pi; % ilk değeri Kappa

XL(ii+1)=XL(ii)+BI(4); % ilk değeri TX

YL(ii+1)=YL(ii)+BI(5); % İlk değeri TYF

ZL(ii+1)=ZL(ii)+BI(6);

end

%% Varyanslar
 % Net Cevaplar
S0=sqrt((V'*V)/(size(B,1)-size(B,2)));
Som=S0*sqrt(Oxx(1,1))*180/pi*3600; % saniye
Sfi=S0*sqrt(Oxx(2,2))*180/pi*3600; % saniye
Ska=S0*sqrt(Oxx(3,3))*180/pi*3600; % saniye
SXL=S0*sqrt(Oxx(4,4))*1000; % mm
SYL=S0*sqrt(Oxx(5,5))*1000; % mm 
SZL=S0*sqrt(Oxx(6,6))*1000; % mm

fprintf("Genel Hata: %f\n",S0);
fprintf("Omega: %f Omega hatası: %f\n",om(ii+1),Som);
fprintf("Phi: %f Phi hatası: %f\n",fi(ii+1),Sfi);
fprintf("Kappa: %f Kappa hatası: %f\n",kap(ii+1),Ska);
fprintf("XL: %f XL hatası: %f\n",XL(ii+1),SXL);
fprintf("YL: %f YL hatası: %f\n",YL(ii+1),SYL);
fprintf("ZL: %f ZL hatası: %f\n",ZL(ii+1),SZL);
