% Solves Ax = b with 
% 
%   A = H+S 
%
% using H^{-1} as the left-preconditioner.
% Hsolve is intented to be a function which implements some H solve, e.g. @(r) H \ r. This is done so one can use inexact Hsolves if desired.
% This implementation corresponds directly to Algorithm 2 in reference [12] of the paper. 
% 
%
% Logan Reed Oct 17, 2026.

function [x, iter] = rapoport(A, H, S, b, Hsolve, tol, maxit)

  a = zeros(1, maxit + 3);

  % Initialization
  x  = zeros(numel(b), 1);
  tv = Hsolve(b);
  a(1) = sqrt(abs(real(b' * tv)));           
  V2 = tv / a(1);

  % Givens state
  g = struct('c', 1, 's', 0, 'cc', 1, 'eps', a(1));   

  iter = 0;
  % the relative residual of x_0 is 1
  if 1 <= tol || maxit < 1
    return;
  end

  % Base steps needed before starting the recurrence
  tv = -Hsolve(S * V2);
  a(2) = sqrt(abs(real(tv' * (H * tv))));
  Vkm1 = tv / a(2);

  tv = -Hsolve(S * Vkm1) + a(2) * V2;
  a(3) = sqrt(abs(real(tv' * (H * tv))));
  Vk = tv / a(3);

  Vx = V2   / (1 + a(2)^2);
  Vz = Vkm1 / (1 + a(2)^2 + a(3)^2);

  x = a(1) * Vx;
  g = givens(g, 0, a(2));
  iter = 1;
  if g.eps / a(1) <= tol || maxit < 2
    return;
  end

  x = a(1) * (Vx + a(2) * Vz);
  g = givens(g, a(2), a(3));
  iter = 2;
  if g.eps / a(1) <= tol || maxit < 3
    return;
  end

  ex  = 1 / (1 + a(2)^2);
  ez  = 1 / (1 + a(2)^2 + a(3)^2);
  ey1 = a(1) / (1 + a(2)^2);
  ey2 = (a(1) * a(2)) / (1 + a(2)^2 + a(3)^2);

  % Main Recurrence
  k = 4;
  while true

    tv   = -Hsolve(S * Vk) + a(k-1) * Vkm1;
    a(k) = sqrt(abs(real(tv' * (H * tv))));
    Vk1  = tv / a(k);

    w = 1 / (1 + a(k)^2 + a(k-1)^2 - ex * (a(k-1)^2) * (a(k-2)^2));

    Vxnew = Vz;
    Vznew = w * (Vk + (a(k-2) * a(k-1)) * Vx);

    x = x + (a(k-1) * a(k-2) * ey1) * Vznew;

    exnew  = ez;
    eznew  = w;
    ey1new = ey2;
    ey2new = a(k-1) * a(k-2) * ey1 * eznew;

    g = givens(g, a(k-1), a(k));
    if g.eps / a(1) <= tol || k - 1 >= maxit
      break;
    end

    Vx = Vxnew;  Vz = Vznew;
    ex = exnew;  ez = eznew;  ey1 = ey1new;  ey2 = ey2new;
    Vkm1 = Vk;   Vk = Vk1;
    k = k + 1;
  end

  iter = k - 1;
end

% givens function for orthogonalization part of the main function
function g = givens(g, aj, aj1)
  m1 = g.cc * aj;                    
  m2 = -g.s * m1 + g.c;            
  r  = hypot(m2, aj1);
  c  = m2 / r;
  s  = -aj1 / r;
  g.eps = abs(s * g.eps);
  g.cc = g.c;
  g.c  = c;    g.s  = s;
end
