% Hava fotoğrafı pozlama noktası koordinatları
clc;clear all;
Xt = 9274.2; % m
Yt = 8292.0; % m
Zt = 1500.1; % m

% Rektifikasyon düzlemi yüksekliği
Zh = 110.0; % m

% Kontrol noktalarının koordinatları [X, Y, Z]
kontrol_noktalari = [
    9300.0, 8300.0, 120.0;
    9250.0, 8280.0, 100.0;
    9350.0, 8320.0, 130.0;
    % Daha fazla nokta eklenebilir
];

% Sonuçları saklamak için matris
sonuc = zeros(size(kontrol_noktalari, 1), 2);

% Hesaplama
for i = 1:size(kontrol_noktalari, 1)
    % Kontrol noktasının koordinatları
    X = kontrol_noktalari(i, 1);
    Y = kontrol_noktalari(i, 2);
    Z = kontrol_noktalari(i, 3);
    
    % Orijinal görüntü koordinatları
    x = X - Xt;
    y = Y - Yt;
    
    % Rektifiye edilmiş koordinatlar
    x_rekt = x * (Zh - Zt) / (Z - Zt);
    y_rekt = y * (Zh - Zt) / (Z - Zt);
    
    % Sonuçlara kaydet
    sonuc(i, :) = [x_rekt, y_rekt];
end

% Sonuçları göster
disp('Rektifiye edilmiş koordinatlar (x'', y''):');
disp(sonuc);
