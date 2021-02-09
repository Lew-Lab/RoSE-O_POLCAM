% 190404 Tianben Ding
% Calculate Fisher information matric of second moments
% Using fully pixelated images

function FIM = calFIMSecondM(B, I, s, backg)

FIM = zeros(6, 6);

FIM(1, 1) = 0.5 .* (s^2) .* (sum(sum((B.aa) ./ (I + backg))));
FIM(2, 2) = 0.5 .* (s^2) .* (sum(sum((B.bb) ./ (I + backg))));
FIM(3, 3) = 0.5 .* (s^2) .* (sum(sum((B.cc) ./ (I + backg))));
FIM(4, 4) = 0.5 .* (s^2) .* (sum(sum((B.dd) ./ (I + backg))));
FIM(5, 5) = 0.5 .* (s^2) .* (sum(sum((B.ee) ./ (I + backg))));
FIM(6, 6) = 0.5 .* (s^2) .* (sum(sum((B.ff) ./ (I + backg))));

FIM(1, 2) = (s^2) .* (sum(sum((B.ab) ./ (I + backg))));
FIM(1, 3) = (s^2) .* (sum(sum((B.ac) ./ (I + backg))));
FIM(1, 4) = (s^2) .* (sum(sum((B.ad) ./ (I + backg))));
FIM(1, 5) = (s^2) .* (sum(sum((B.ae) ./ (I + backg))));
FIM(1, 6) = (s^2) .* (sum(sum((B.af) ./ (I + backg))));

FIM(2, 3) = (s^2) .* (sum(sum((B.bc) ./ (I + backg))));
FIM(2, 4) = (s^2) .* (sum(sum((B.bd) ./ (I + backg))));
FIM(2, 5) = (s^2) .* (sum(sum((B.be) ./ (I + backg))));
FIM(2, 6) = (s^2) .* (sum(sum((B.bf) ./ (I + backg))));

FIM(3, 4) = (s^2) .* (sum(sum((B.cd) ./ (I + backg))));
FIM(3, 5) = (s^2) .* (sum(sum((B.ce) ./ (I + backg))));
FIM(3, 6) = (s^2) .* (sum(sum((B.cf) ./ (I + backg))));

FIM(4, 5) = (s^2) .* (sum(sum((B.de) ./ (I + backg))));
FIM(4, 6) = (s^2) .* (sum(sum((B.df) ./ (I + backg))));

FIM(5, 6) = (s^2) .* (sum(sum((B.ef) ./ (I + backg))));

FIM = FIM + FIM.';
end
