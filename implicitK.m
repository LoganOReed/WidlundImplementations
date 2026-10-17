function [x, iter] = implicitK(A, H, S, b, Hsolve, tol, maxit)

  n = numel(b);

  x     = zeros(n, 1);
  xprev = zeros(n, 1);
  r     = b;

  aprev = 0;
  pprev = 1;

  for m = 1:(maxit + 1)

    v = Hsolve(r);
    p = real(r' * v);

    if m == 1
      scale = sqrt(abs(p));
    end
    if sqrt(abs(p)) / scale <= tol || m > maxit
      break;
    end

    if ~isfinite(p) || p <= 0
      break;
    end

    a = 1 + (p / pprev) * aprev;
    w = 1 / a;

    xnext = (1 - w) * xprev + w * x + w * v;
    r     = b - A * xnext;

    aprev = a;
    pprev = p;
    xprev = x;
    x     = xnext;
  end

  iter = m - 1;
end
