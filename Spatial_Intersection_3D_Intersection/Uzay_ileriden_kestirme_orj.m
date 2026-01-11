clear all;
delete *.asv;
format long g;
clc;

% Burada yöneltme parametreleri veriliyor, sonuçta yeryüzü koord.ları
% buluyoruz. Öncesinde geriden kestirme yapıp buraya geçilebilir.

f=152.057;
f=f/1000;

photo=xlsread('Koordinatlar_Uzay_İleriden_Kestirme.xlsx',"Örnek 11-2");
OM=photo(:,1); 
FI=photo(:,2); 
KAP=photo(:,3);
XL=photo(:,4);
YL=photo(:,5);
ZL=photo(:,6);

xy=xlsread('Koordinatlar_Uzay_İleriden_Kestirme.xlsx','Örnek 11-2 Noktalar');
x1=xy(:,1).*10^-3;
y1=xy(:,2).*10^-3;
x2=xy(:,3).*10^-3;
y2=xy(:,4).*10^-3;
ng=size(xy,1);

for i=1:3*ng;
    Delta(i,1)=1;
end


B=sqrt((XL(2)-XL(1))^2+(YL(2)-YL(1))^2);
H=mean(ZL)

for i=1:ng;
    pa(i,1)=x1(i,1)-x2(i,1);
    xaussu(i,1)=B*x1(i,1)/pa(i,1);
    yaussu(i,1)=B*y1(i,1)/pa(i,1);
    ZA(i,1)=H-B*f/pa(i,1);
    XA(i,1)=((XL(2)-XL(1))*xaussu(i,1)/B)-((YL(2)-YL(1))*yaussu(i,1)/B)+XL(1);
    YA(i,1)=((XL(2)-XL(1))*yaussu(i,1)/B)+((YL(2)-YL(1))*xaussu(i,1)/B)+YL(1);

    XXXA(i,1)=XA(i,1);
    YYYA(i,1)=YA(i,1);
    ZZZA(i,1)=ZA(i,1);
end
disp(XA);
disp(YA);
disp(ZA);
% m1 matrisi için
m1(1,1) = cosd(FI(1)) * cosd(KAP(1));
m1(1,2) = sind(OM(1)) * sind(FI(1)) * cosd(KAP(1)) + cosd(OM(1)) * sind(KAP(1));
m1(1,3) = -cosd(OM(1)) * sind(FI(1)) * cosd(KAP(1)) + sind(OM(1)) * sind(KAP(1));

m1(2,1) = -cosd(FI(1)) * sind(KAP(1));
m1(2,2) = -sind(OM(1)) * sind(FI(1)) * sind(KAP(1)) + cosd(OM(1)) * cosd(KAP(1));
m1(2,3) = cosd(OM(1)) * sind(FI(1)) * sind(KAP(1)) + sind(OM(1)) * cosd(KAP(1));

m1(3,1) = sind(FI(1));
m1(3,2) = -sind(OM(1)) * cosd(FI(1));
m1(3,3) = cosd(OM(1)) * cosd(FI(1));

% m2 matrisi için
m2(1,1) = cosd(FI(2)) * cosd(KAP(2));
m2(1,2) = sind(OM(2)) * sind(FI(2)) * cosd(KAP(2)) + cosd(OM(2)) * sind(KAP(2));
m2(1,3) = -cosd(OM(2)) * sind(FI(2)) * cosd(KAP(2)) + sind(OM(2)) * sind(KAP(2));

m2(2,1) = -cosd(FI(2)) * sind(KAP(2));
m2(2,2) = -sind(OM(2)) * sind(FI(2)) * sind(KAP(2)) + cosd(OM(2)) * cosd(KAP(2));
m2(2,3) = cosd(OM(2)) * sind(FI(2)) * sind(KAP(2)) + sind(OM(2)) * cosd(KAP(2));

m2(3,1) = sind(FI(2));
m2(3,2) = -sind(OM(2)) * cosd(FI(2));
m2(3,3) = cosd(OM(2)) * cosd(FI(2));

ii=0;
while max(abs(Delta))>0.001
    ii=ii+1;
    disp([num2str(ii),' nci iterasyon']);

a=1;
q=1;
k=3;

for i=1:1:ng
    % r1, s1, q1 için (m1 matrisi kullanılarak)
    r1(i,1) = m1(1,1)*(XA(i,1)-XL(1)) + m1(1,2)*(YA(i,1)-YL(1)) + m1(1,3)*(ZA(i,1)-ZL(1));
    s1(i,1) = m1(2,1)*(XA(i,1)-XL(1)) + m1(2,2)*(YA(i,1)-YL(1)) + m1(2,3)*(ZA(i,1)-ZL(1));
    q1(i,1) = m1(3,1)*(XA(i,1)-XL(1)) + m1(3,2)*(YA(i,1)-YL(1)) + m1(3,3)*(ZA(i,1)-ZL(1));

    % r2, s2, q2 için (m2 matrisi kullanılarak)
    r2(i,1) = m2(1,1)*(XA(i,1)-XL(2)) + m2(1,2)*(YA(i,1)-YL(2)) + m2(1,3)*(ZA(i,1)-ZL(2));
    s2(i,1) = m2(2,1)*(XA(i,1)-XL(2)) + m2(2,2)*(YA(i,1)-YL(2)) + m2(2,3)*(ZA(i,1)-ZL(2));
    q2(i,1) = m2(3,1)*(XA(i,1)-XL(2)) + m2(3,2)*(YA(i,1)-YL(2)) + m2(3,3)*(ZA(i,1)-ZL(2));


    b(a,q)=f/q1(i,1)^2*(r1(i,1)*m1(3,1)-q1(i,1)*m1(1,1));
    b(a,q+1)=f/q1(i,1)^2*(r1(i,1)*m1(3,2)-q1(i,1)*m1(1,2));
    b(a,q+2)=f/q1(i,1)^2*(r1(i,1)*m1(3,3)-q1(i,1)*m1(1,3));

    b(a+1,q)=f/q1(i,1)^2*(s1(i,1)*m1(3,1)-q1(i,1)*m1(2,1));
    b(a+1,q+1)=f/q1(i,1)^2*(s1(i,1)*m1(3,2)-q1(i,1)*m1(2,2));
    b(a+1,q+2)=f/q1(i,1)^2*(s1(i,1)*m1(3,3)-q1(i,1)*m1(2,3));

    b(k,q)=f/q2(i,1)^2*(r2(i,1)*m2(3,1)-q2(i,1)*m2(1,1));
    b(k,q+1)=f/q2(i,1)^2*(r2(i,1)*m2(3,2)-q2(i,1)*m2(1,2));
    b(k,q+2)=f/q2(i,1)^2*(r2(i,1)*m2(3,3)-q2(i,1)*m2(1,3));

    b(k+1,q)=f/q2(i,1)^2*(s2(i,1)*m2(3,1)-q2(i,1)*m2(2,1));
    b(k+1,q+1)=f/q2(i,1)^2*(s2(i,1)*m2(3,2)-q2(i,1)*m2(2,2));
    b(k+1,q+2)=f/q2(i,1)^2*(s2(i,1)*m2(3,3)-q2(i,1)*m2(2,3));

    J1(i,1)=x1(i,1)+f*r1(i,1)/q1(i,1); 
    K1(i,1)=y1(i,1)+f*s1(i,1)/q1(i,1);

    J2(i,1)=x2(i,1)+f*r2(i,1)/q2(i,1);
    K2(i,1)=x2(i,1)+f*s2(i,1)/q2(i,1);

    E(a,1)=J1(i,1);
    E(a+1,1)=K1(i,1);

    E(k,1)=J2(i,1);
    E(k+1,1)=K2(i,1);

    a=a+4;
    k=k+4;
    q=q+3;
end


Qxx=inv(b'*b);
BI=Qxx*(b'*E);
V=b*BI-E;

Delta=BI(:,1);

t=1;

for i=1:ng;
    XA(i,1)=XA(i,1)+BI(t);
    YA(i,1)=YA(i,1)+BI(t+1);
    ZA(i,1)=ZA(i,1)+BI(t+2);

t=t+3;
end 

XXXA(:,ii+1)=XA(:,1);
YYYA(:,ii+1)=YA(:,1);
ZZZA(:,ii+1)=ZA(:,1);

end


SO=sqrt((V'*V)/(size(B,1)-size(B,2)));


[num,txt,raw]=xlsread('Koordinatlar_Uzay_İleriden_Kestirme.xlsx','Örnek 11-2 Noktalar');
a=1;
for i=1:ng;
    T(a)=char(txt(i+2,1));
    a=a+1;
end

fileID=fopen("Uzay_ileriden_kestirme_Örnek 11-2.txt","w");

fprintf(fileID,'%16s\n \n','Nesne Uzay Koordinatları:');
fprintf(fileID,'%0s %11s %16s %17s\n','Nokta','X,mm','Y,mm','Z,mm');

for i=1:ng;
    fprintf(fileID,'%c %15.4f %17.4f %17.4f\n',T(i),XXXA(i,ii+1),YYYA(i,ii+1),ZZZA(i,ii+1));
end

fprintf(fileID,'\n');
% fclose(fileID);


















