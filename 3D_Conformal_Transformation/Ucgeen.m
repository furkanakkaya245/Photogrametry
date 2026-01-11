function [hmax, satir, sutun, uX, uY, uZ, ux, uy, uz] = Ucgeen(Koordinatlar)
    [x, y] = size(Koordinatlar);

    % Koordinatlar matrisinin yeterli satıra sahip olup olmadığını kontrol et
    if x < 3
        error('Koordinatlar matrisinde en az 3 satır olmalıdır.');
    end
    
    % 3'lü kombinasyonları bulma
    uX = nchoosek(Koordinatlar(:, 1), 3);
    uY = nchoosek(Koordinatlar(:, 2), 3);
    uZ = nchoosek(Koordinatlar(:, 3), 3);
    ux = nchoosek(Koordinatlar(:, 4), 3);
    uy = nchoosek(Koordinatlar(:, 5), 3);
    uz = nchoosek(Koordinatlar(:, 6), 3);

    % Kenar uzunluklarını hesapla
    AB = zeros(size(uX, 1), 3);  % Kenarları saklayacak matrisin boyutunu uX'in satır sayısına göre ayarla
    for i = 1:size(uX, 1)  % uX'teki kombinasyon sayısı kadar döngü
        AB(i, 1) = sqrt((uX(i, 2) - uX(i, 1))^2 + (uY(i, 2) - uY(i, 1))^2);
        AB(i, 2) = sqrt((uX(i, 3) - uX(i, 1))^2 + (uY(i, 3) - uY(i, 1))^2);
        AB(i, 3) = sqrt((uX(i, 3) - uX(i, 2))^2 + (uY(i, 3) - uY(i, 2))^2);     
    end

    % Her üçgen için yarı çevreyi hesapla
    toplam = sum(AB, 2);
    s = toplam / 2;

    % Her üçgenin alanını hesapla
    Alan = zeros(length(s), 1);
    for i = 1:length(s)
        Alan(i) = sqrt(s(i) * (s(i) - AB(i, 1)) * (s(i) - AB(i, 2)) * (s(i) - AB(i, 3)));
    end

    % Her üçgenin yüksekliğini hesapla
    h = zeros(size(uX, 1), 1);
    for i = 1:size(uX, 1)
        h(i) = 2 * Alan(i) / max(AB(i, 1:3));  % Her üçgenin yüksekliğini en uzun kenara göre hesapla
    end

    % En yüksek yüksekliğin bulunması
    [hmax, idx] = max(h);  % max(h) ile en yüksek yüksekliği bul
    [satir, sutun] = ind2sub(size(uX), idx);  % İndeksi satır ve sütun formatına dönüştür
end
