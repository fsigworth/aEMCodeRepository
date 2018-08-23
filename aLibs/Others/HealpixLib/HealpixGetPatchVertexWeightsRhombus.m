% Get interpolation weights for the rhombus patch
% vertexes : v1(0, -1), v2(0.5, 0), v3(0, 1), v4(-0.5, 0)
% [w1, w2, w3, w4] = HealpixGetPatchVertexWeightsRhombus(delta_i, delta_j)
%
% Parameters
% delta_i, delta_j : position in rhombus (center is (0, 0))
% w1 : weight for v1 (north)
% w2 : weight for v3 (south
% w3 : weight for v4 (west)
% w4 : weight for v2 (east)

function [w1, w2, w3, w4] = HealpixGetPatchVertexWeightsRhombus(delta_i, delta_j)

% Let v1(0, -1), v2(0.5, 0), v3(0, 1), v4(-0.5, 0)
sn = [delta_i, delta_j] * [ 2; 1] / sqrt(5);	% v1->v2 �x�N�g������
so = [delta_i, delta_j] * [-1; 2] / sqrt(5);	% v1->v2 �x�N�g���ɒ��s�ȃx�N�g������

h0 = sqrt(5) / 5 - so;  % v1���܂ޕ��s�l�ӌ`��v1-v2���ӂƂ����ۂ̍���
h1 = sqrt(5) / 5 + so;  % v4���܂ޕ��s�l�ӌ`��v3-v4���ӂƂ����ۂ̍���

% (delta_i, delta_j)��ʂ�v2-v3�ɕ��s�Ȓ�����v1-v2�̌�_�����߂�
px = (2 * delta_j + delta_i + 1) / 4;
py = -2 * px - 1;
% (px, py)��v1�̒���
q = sqrt(px^2 + (py + 1)^2);

b0 = q;                 % v1���܂ޕ��s�l�ӌ`�̒��(v1-v2���ɉ�������)
b1 = sqrt(5/4) - b0;    % v2���܂ޕ��s�l�ӌ`�̒��(v1-v2���ɉ�������)

% north
w1 = b1 * h1;   % v3���܂ޕ��s�l�ӌ`�̖ʐ�

% south
w2 = b0 * h0;   % v1���܂ޕ��s�l�ӌ`�̖ʐ�

% west
w3 = b1 * h0;   % v2���܂ޕ��s�l�ӌ`�̖ʐ�

% east
w4 = b0 * h1;   % v4���܂ޕ��s�l�ӌ`�̖ʐ�
