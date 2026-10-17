function [x, iter] = mss(A, H, S, b, Hsolve, tol, maxit)

  n = numel(b);

  bhat  = Hsolve(b);
  beta0 = sqrt(abs(real(b' * bhat)));

  v     = bhat / beta0;
  vprev = zeros(n, 1);
  betam = 0;
  p     = zeros(n, 1);

  x   = zeros(n, 1);
  est = 1;

  for m = 1:(maxit + 1)

    if est <= tol || m > maxit
      break;
    end

    w = Hsolve(S * v);
    if m > 1
      w = w + betam * vprev;
    end

    beta = sqrt(abs(real(w' * (H * w))));

    if m == 1
      eta  = 1;
      zeta = beta0;
    else
      lam  = betam / eta;
      eta  = 1 + lam * betam;
      zeta = -lam * zeta;
    end
    p = (v + betam * p) / eta;
    x = x + zeta * p;

    est = beta * abs(zeta / eta) / beta0;

    if beta <= eps * beta0
      iter = m;
      return;
    end

    vprev = v;
    v     = w / beta;
    betam = beta;
  end

  iter = m - 1;
end

