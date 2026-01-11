clc ; clear all ; close all; format long g ; delete *.asv;
% Gercek yeryuzu koordinatlarını bulun..
%% Dosyadan Okuma

% Bu kısımda en yüksek irtifaya sahip ücgeni olusturan noktaları ve
% bunun yuksekligini bulduk.

Koordinatlar = xlsread('mutlak_yoneltme.xlsx')
[x,y] = size(Koordinatlar);

[hmax,satir,sutun,uX,uY,uZ,ux,uy,uz]=Ucgeen(Koordinatlar);



for i=1:3;
   
    j=1:x;
    
    XYZ(i,j)=Koordinatlar(j,i);

    xyz(i,j)=Koordinatlar(j,i+3);

end
delta1 = [1 1 1 ];delta2 = [1 1 1 ];
%% Normal Denklemler
% Keyfi koordinatlar için normal vektor Na
% Kontrol koordinatları için normal vektor Nc(yeryuzu nokt.)

Na=[(uy(satir,2)-uy(satir,1))*(uz(satir,3)-uz(satir,1))-(uy(satir,3)-uy(satir,1))*(uz(satir,2)-uz(satir,1))
    (ux(satir,3)-ux(satir,1))*(uz(satir,2)-uz(satir,1))-(ux(satir,2)-ux(satir,1))*(uz(satir,3)-uz(satir,1))
    (ux(satir,2)-ux(satir,1))*(uy(satir,3)-uy(satir,1))-(ux(satir,3)-ux(satir,1))*(uy(satir,2)-uy(satir,1))];
 

Nc=[(uY(satir,2)-uY(satir,1))*(uZ(satir,3)-uZ(satir,1))-(uY(satir,3)-uY(satir,1))*(uZ(satir,2)-uZ(satir,1))
    (uX(satir,3)-uX(satir,1))*(uZ(satir,2)-uZ(satir,1))-(uX(satir,2)-uX(satir,1))*(uZ(satir,3)-uZ(satir,1))
    (uX(satir,2)-uX(satir,1))*(uY(satir,3)-uY(satir,1))-(uX(satir,3)-uX(satir,1))*(uY(satir,2)-uY(satir,1))];

% Burada her normal vektör için eğim ve azimut hesaplanıyor.(derece)

ta=atand(Na(3)/sqrt(Na(1)^2+Na(2)^2))+90    %tilt_a
aa=atand(Na(1)/Na(2))-180                   %azimuth_a
tc=atand(Nc(3)/sqrt(Nc(1)^2+Nc(2)^2))+90    %tilt_c
ac=atand(Nc(1)/Nc(2))-180                   %azimuth_c
sa=0;
sc=0;
% Salınım sıfır olarak ayarlanmıştır..
% Aşağıda iki tane döndürme matrisi oluşturuluyor.
% Matrisler 101 ve 102 noktalarına uyarlanarak x-y düzleminde yatay.. 
% iki nokta çiftinin koord.ları elde ediliyor.

ma11=-cosd(aa)*cosd(sa)-sind(aa)*cosd(ta)*sind(sa);
ma12=sind(aa)*cosd(sa)-cosd(aa)*cosd(ta)*sind(sa);
ma13=-sind(ta)*sind(sa);
ma21=cosd(aa)*sind(sa)-sind(aa)*cosd(ta)*cosd(sa);
ma22=-sind(aa)*sind(sa)-cosd(aa)*cosd(ta)*cosd(sa);
ma23=-sind(ta)*cosd(sa);
ma31=-sind(aa)*sind(ta);
ma32=-cosd(aa)*sind(ta);
ma33=cosd(ta);

Ma1=[ma11 ma12 ma13;ma21 ma22 ma23;ma31 ma32 ma33]


mc11=-cosd(ac)*cosd(sc)-sind(ac)*cosd(tc)*sind(sc);
mc12=sind(ac)*cosd(sc)-cosd(ac)*cosd(tc)*sind(sc);
mc13=-sind(tc)*sind(sc);
mc21=cosd(ac)*sind(sc)-sind(ac)*cosd(tc)*cosd(sc);
mc22=-sind(ac)*sind(sc)-cosd(ac)*cosd(tc)*cosd(sc);
mc23=-sind(tc)*cosd(sc);
mc31=-sind(ac)*sind(tc);
mc32=-cosd(ac)*sind(tc);
mc33=cosd(tc);

Mc=[mc11 mc12 mc13;mc21 mc22 mc23;mc31 mc32 mc33]

% Donukluk matrisini noktaların koordinatları ile carpiyoruz ve dondurulmus
% noktaları elde ediyoruz.
% Ka --> keyfi sistemde dondurulmus noktalar
% Kc --> kontrol sisteminde dondurulmus noktalar
Ka101 = Ma1*[ux(satir,1);uy(satir,1);uz(satir,1)]
Ka102 = Ma1*[ux(satir,2);uy(satir,2);uz(satir,2)]
 
aza=atan2((Ka102(1)-Ka101(1)),(Ka102(2)-Ka101(2)))*180/pi   %azimuth_arbitrary(gecici)

Kc101 = Mc*[uX(satir,1);uY(satir,1);uZ(satir,1)]
Kc102 = Mc*[uX(satir,2);uY(satir,2);uZ(satir,2)]

azc=atan2((Kc102(1)-Kc101(1)),(Kc102(2)-Kc101(2)))*180/pi   %azimuth_control

% Dondurulmus cizgiler arasındaki salınım(sa), dondurulen noktalar
% arasındaki azimut farkları kullanılarak hesaplanır.
sa=azc-aza            %swing_a

%% Donukluk Matrisi

mb11=-cosd(aa)*cosd(sa)-sind(aa)*cosd(ta)*sind(sa);
mb12=sind(aa)*cosd(sa)-cosd(aa)*cosd(ta)*sind(sa);
mb13=-sind(ta)*sind(sa);
mb21=cosd(aa)*sind(sa)-sind(aa)*cosd(ta)*cosd(sa);
mb22=-sind(aa)*sind(sa)-cosd(aa)*cosd(ta)*cosd(sa);
mb23=-sind(ta)*cosd(sa);
mb31=-sind(aa)*sind(ta);
mb32=-cosd(aa)*sind(ta);
mb33=cosd(ta);

Ma2=[mb11 mb12 mb13;mb21 mb22 mb23;mb31 mb32 mb33]

% Mac genel donukluk matrisidir. Yani iki sistem arasındaki donukluk matrisi.
Mac=Ma2'*Mc

% Omega phi kappa değerlerinin elde edilmesi
% bunlar ilk yaklaşık değerler olarak geçer
OM=wrapTo2Pi(atan2(-Mac(3,2),Mac(3,3)))*180/pi 
FI=asind(Mac(3,1))
KA=wrapTo2Pi(atan2(-Mac(2,1),Mac(1,1)))*180/pi 


%% Yaklasik Olcegin Hesaplanmasi (s)
a=1;
for i=1:1:x;
    for j=1:1:x;
        if i==j || j<i;
            continue;
        else
            AB (a)=sqrt((Koordinatlar(j,1)-Koordinatlar(i,1))^2+(Koordinatlar(j,2)-Koordinatlar(i,2))^2);
            ab (a)=sqrt((Koordinatlar(j,4)-Koordinatlar(i,4))^2+(Koordinatlar(j,5)-Koordinatlar(i,5))^2);
            a=a+1;
        end
    end
end

for i=1:length(AB);
    ss(i)=AB(i)/ab(i);
end

s=mean(ss)

%% Otelemelerin Hesaplanmasi

M11=cosd(FI)*cosd(KA);
M12=sind(OM)*sind(FI)*cosd(KA)+cosd(OM)*sind(KA);
M13=-cosd(OM)*sind(FI)*cosd(KA)+sind(OM)*sind(KA);
M21=-cosd(FI)*sind(KA);
M22=-sind(OM)*sind(FI)*sind(KA)+cosd(OM)*cosd(KA);
M23=cosd(OM)*sind(FI)*sind(KA)+sind(OM)*cosd(KA);
M31=sind(FI);
M32=-sind(OM)*cosd(FI);
M33=cosd(OM)*cosd(FI);

M = [M11, M12, M13; M21, M22, M23; M31, M32, M33];

Txyz=XYZ-s*M'*xyz;

Tx=mean(Txyz(1,:))
Ty=mean(Txyz(2,:))
Tz=mean(Txyz(3,:))



 ii=0;
 while max(abs(delta1)) >0.0000005 || max(abs(delta2)) > 0.0001
    ii=ii+1;    
     
      
disp([num2str(ii),' nci iterasyon.'])


%% Katsayılar Matrisi
 
M11=cosd(FI)*cosd(KA);
M12=sind(OM)*sind(FI)*cosd(KA)+cosd(OM)*sind(KA);
M13=-cosd(OM)*sind(FI)*cosd(KA)+sind(OM)*sind(KA);
M21=-cosd(FI)*sind(KA);
M22=-sind(OM)*sind(FI)*sind(KA)+cosd(OM)*cosd(KA);
M23=cosd(OM)*sind(FI)*sind(KA)+sind(OM)*cosd(KA);
M31=sind(FI);
M32=-sind(OM)*cosd(FI);
M33=cosd(OM)*cosd(FI);


a=1;
for i=1:x;
    b(a,1)=M11*xyz(1,i)+M21*xyz(2,i)+M31*xyz(3,i);
    b(a,2)=0;
    b(a,3)=((-sind(FI)*cosd(KA))*xyz(1,i)+(sind(FI)*sind(KA))*xyz(2,i)+(cosd(FI))*xyz(3,i))*s;
    b(a,4)=(M21*xyz(1,i)-M11*xyz(2,i))*s;
    b(a,5)=1;
    b(a,6)=0;
    b(a,7)=0;
    b(a+1,1)=M12*xyz(1,i)+M22*xyz(2,i)+M32*xyz(3,i);
    b(a+1,2)=(-M13*xyz(1,i)-M23*xyz(2,i)-M33*xyz(3,i))*s;   
    b(a+1,3)=((sind(OM)*M11)*xyz(1,i)+(sind(OM)*M21)*xyz(2,i)+(sind(OM)*M31)*xyz(3,i))*s;
    b(a+1,4)=(M(2,2)*xyz(1,i)-M(1,2)*xyz(2,i))*s;
    b(a+1,5)=0;
    b(a+1,6)=1;
    b(a+1,7)=0;
    b(a+2,1)=M13*xyz(1,i)+M23*xyz(2,i)+M33*xyz(3,i);
    b(a+2,2)=(M12*xyz(1,i)+M22*xyz(2,i)+M32*xyz(3,i))*s;
    b(a+2,3)=((-cosd(OM)*M11)*xyz(1,i)+(-cosd(OM)*M21)*xyz(2,i)+(-cosd(OM)*M31)*xyz(3,i))*s;
    b(a+2,4)=(M23*xyz(1,i)-M13*xyz(2,i))*s;
    b(a+2,5)=0;
    b(a+2,6)=0;
    b(a+2,7)=1;
     
    a=a+3;
end

A=[b(:,1) b(:,2) b(:,3) b(:,4) b(:,5) b(:,6) b(:,7)];
    
 %% L Matrisi
 a=1;
for i=1:x;
    L(a,1)=[XYZ(1,i)]-Tx;
    L(a+1,1)=[XYZ(2,i)]-Ty;
    L(a+2,1)=[XYZ(3,i)]-Tz;
    a=a+3;
end
     
%% Dengeleme

BI = inv(A'*A)*A'*L
V = A*BI-L

s1=BI(1,1)

delta1 = BI(2:4)
delta2 = BI(5:7)

OM = OM+(delta1(1,1))*180/pi
FI = FI+(delta1(2,1))*180/pi
KA = KA+(delta1(3,1))*180/pi
Tx = Tx+delta2(1,1)
Ty = Ty+delta2(2,1)
Tz = Tz+delta2(3,1)



 end 

 disp('******************Sonuçlar******************');
 disp('*****BI*****')
 disp(BI);
 disp('*****V*****')
 disp(V);
 disp('*****L*****')
 disp(L);
 disp('*****s1*****')
 disp(s1);
 disp('**********Delta-1**********');
 disp(delta1);
 disp('**********Delta-2**********');
 disp(delta2);
 
 


 disp('DÖNÜŞÜM PARAMETRELERİ: ');
 disp(['Ölçek Katsayısı(°) ',num2str(s)]);
 disp(['Dengelenmis Omega(°) ',num2str(OM)]);
 disp(['Dengelenmis Phi(°) ',num2str(FI)]);
 disp(['Dengelenmis Kappa(°) ',num2str(KA)]);
 disp(['Dengelenmiş Tx(m) ',num2str(Tx)]);
 disp(['Dengelenmiş Ty(m) ',num2str(Ty)]);
 disp(['Dengelenmiş Tz(m) ',num2str(Tz)]);



 corA = [[-4.8352,1.9730,1.0888],
         [89.0970,2.7047,0.3391],
         [89.2672,82.8667,1.7862],
         [0,0,152.113],
         [91.9740,-1.7346,148.3015]
         ];


 sonuc = [[9265.105,10213.339,64.073],
          [9575.295,10220.215,66.213],
          [9572.011,10485.010,66.406],
          [9273.552,10215.603,563.122],
          [9577.546,10214.067,555.197];
 ];
 
 sonuc./corA;

delta = [OM, FI, KA, Tx, Ty, Tz,s]; % Delta ve A aynı boyutta olmalı
F = delta .* A; % Eleman bazlı çarpım
result = inv( A' * A) * F' % F transpozu alınarak matris çarpımı yapılabilir


Msonuc=hesapla_M(OM,FI,KA);
Noktalar = [-4.8352, 89.0970, 89.2672, 0, 91.9740;
             1.9730,  2.7047, 82.8667, 0, -1.7346;
             1.0888,  0.3391,  1.7862, 152.113, 148.3015];
T = [Tx; Ty; Tz]; 
YerKoor = s * (Msonuc' * Noktalar) + T;
YerKoor_Liste = YerKoor'







function M = hesapla_M(omega_deg, phi_deg, kappa_deg)
    w = deg2rad(omega_deg);
    p = deg2rad(phi_deg);
    k = deg2rad(kappa_deg);
    sw = sin(w); cw = cos(w);
    sp = sin(p); cp = cos(p);
    sk = sin(k); ck = cos(k);
    M = zeros(3,3);
    
    M(1,1) = cp * ck;
    M(1,2) = sw * sp * ck + cw * sk;
    M(1,3) = -cw * sp * ck + sw * sk;
    
    M(2,1) = -cp * sk;
    M(2,2) = -sw * sp * sk + cw * ck;
    M(2,3) = cw * sp * sk + sw * ck;
    
    M(3,1) = sp;
    M(3,2) = -sw * cp;
    M(3,3) = cw * cp;
end



