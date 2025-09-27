function u=esc(t,To)
%
%       u=esc(t,To)
% Genera la señal escalón unitario en tiempo discreto
% con retardo mayor o igual a cero.
%
%    u  - señal escalón.
%    t  - vector de tiempo discreto.
%    To - retardo en tiempo continuo, (To>=0).
%    No - retardo en tiempo discreto, No>=0.
%
[Dum N]=size(t);
No=ceil(abs(To)/t(2));
u1= zeros([1 No]);
u2= ones([1 N-No]);
u=[u1 u2];