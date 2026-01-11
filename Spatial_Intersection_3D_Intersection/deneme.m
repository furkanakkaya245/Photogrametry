clc; clear; close all; format long g;

%% 1. Verilenler ve Kamera Parametreleri
f = 152.057; % mm
cam1.pos = [9577.252, 10214.285, 555.192]; % XL1, YL1, ZL1
cam1.angles = [1.4022, -0.3112, 0.6470];   % omega, phi, kappa (derece)

cam2.pos = [9803.241, 10219.622, 556.601]; % XL2, YL2, ZL2
cam2.angles = [-0.1557, -1.7063, 0.5513];  % omega, phi, kappa (derece)

% Fotoğraf Koordinatları (mm) [x1, y1, x2, y2]
noktalar = {
    'A', [-12.843, -54.155, -87.550, -51.157]
    'B', [-22.720, -4.828,  -98.272, -1.702]
    'C', [-0.433,  70.321,  -75.465, 75.310]
    'D', [-25.993, -58.086, -101.077, -55.213]
};

fprintf('%-8s %-15s %-15s %-15s\n', 'Nokta', 'X (m)', 'Y (m)', 'Z (m)');
disp('------------------------------------------------------------');

%% 2. Her Nokta İçin Hesaplama Döngüsü
for i = 1:size(noktalar, 1)
    isim = noktalar{i,1};
    coords = noktalar{i,2};
    
    x1 = coords(1)/1000; y1 = coords(2)/1000;
    x2 = coords(3)/1000; y2 = coords(4)/1000;
    f_m = f/1000;
    
    % --- Başlangıç Yaklaşık Değerleri ---
    % Paralaks mantığıyla kabaca bir yer tahmini yapıyoruz (Diverjansı önlemek için)
    X_est = (cam1.pos(1) + cam2.pos(1)) / 2;
    Y_est = (cam1.pos(2) + cam2.pos(2)) / 2;
    Z_est = 65.0; % Bölüm 11'deki örnek arazi yüksekliği
    
    % Rotasyon Matrisleri
    M1 = hesapla_M(cam1.angles);
    M2 = hesapla_M(cam2.angles);
    
    % --- İterasyon (Newton-Raphson) ---
    for iter = 1:15
        A = []; L = [];
        
        % Foto 1 ve Foto 2 için denklemleri kur
        % A = [B14 B15 B16] matrisi (Kısmi Türevler)
        [A1, L1] = kolinearite_satiri(f_m, M1, cam1.pos, [x1, y1], [X_est, Y_est, Z_est]);
        [A2, L2] = kolinearite_satiri(f_m, M2, cam2.pos, [x2, y2], [X_est, Y_est, Z_est]);
        
        A = [A1; A2];
        L = [L1; L2];
        
        % En Küçük Kareler Çözümü: delta = (A'A) \ (A'L)
        dX = (A' * A) \ (A' * L);
        
        X_est = X_est + dX(1);
        Y_est = Y_est + dX(2);
        Z_est = Z_est + dX(3);
        
        if max(abs(dX)) < 0.0001, break; end
    end
    
    fprintf('%-8s %-15.4f %-15.4f %-15.4f\n', isim, X_est, Y_est, Z_est);
end

%% Yardımcı Fonksiyonlar
function M = hesapla_M(angles)
    w = deg2rad(angles(1)); p = deg2rad(angles(2)); k = deg2rad(angles(3));
    M(1,1) = cos(p)*cos(k);
    M(1,2) = sin(w)*sin(p)*cos(k) + cos(w)*sin(k);
    M(1,3) = -cos(w)*sin(p)*cos(k) + sin(w)*sin(k);
    M(2,1) = -cos(p)*sin(k);
    M(2,2) = -sin(w)*sin(p)*sin(k) + cos(w)*cos(k);
    M(2,3) = cos(w)*sin(p)*sin(k) + sin(w)*sin(k);
    M(3,1) = sin(p);
    M(3,2) = -sin(w)*cos(p);
    M(3,3) = cos(w)*cos(p);
end

function [A, L] = kolinearite_satiri(f, M, Lpos, photo, ground)
    dX = ground(1) - Lpos(1);
    dY = ground(2) - Lpos(2);
    dZ = ground(3) - Lpos(3);
    
    r = M(1,1)*dX + M(1,2)*dY + M(1,3)*dZ;
    s = M(2,1)*dX + M(2,2)*dY + M(2,3)*dZ;
    q = M(3,1)*dX + M(3,2)*dY + M(3,3)*dZ;
    
    % Jacobian Katsayıları (B14, B15, B16)
    A(1,1) = (f/q^2)*(r*M(3,1) - q*M(1,1));
    A(1,2) = (f/q^2)*(r*M(3,2) - q*M(1,2));
    A(1,3) = (f/q^2)*(r*M(3,3) - q*M(1,3));
    
    A(2,1) = (f/q^2)*(s*M(3,1) - q*M(2,1));
    A(2,2) = (f/q^2)*(s*M(3,2) - q*M(2,2));
    A(2,3) = (f/q^2)*(s*M(3,3) - q*M(2,3));
    
    % L Vektörü (Gözlem - Hesaplanan)
    L(1,1) = -(photo(1) + f * r / q);
    L(2,1) = -(photo(2) + f * s / q);
end