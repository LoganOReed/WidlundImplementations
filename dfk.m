function [x, iter] = dfk(A, H, S, b, Hsolve, tol, maxit)

  n = numel(b);

  x = zeros(n, 1);
  r = b;                             

  v_old = zeros(n, 1);               
  z_old = zeros(n, 1);               

  z     = Hsolve(r);                 
  beta0 = sqrt(abs(z' * r));         

  v = (1 / beta0) * r;               
  z = (1 / beta0) * z;               

  y_old     = zeros(n, 1);           
  y_veryold = zeros(n, 1);           

  % the two most recent Givens rotations
  s_veryold = 0;  c_veryold = 1;
  s_old     = 0;  c_old     = 1;

  % trailing entry of Q'(beta0 e_1), so far Q = I
  epsilont     = beta0;
  epsilon      = 0;
  epsilont_old = 0;

  alphat_old  = 0;
  alphatt_old = 1;

  % ||r_0||_{H^-1} / ||b||_{H^-1}
  est = 1;

  for m = 1:(maxit + 1)

    if est <= tol || m > maxit
      break;
    end

    % Next Lanczos step
    w     = A * z;
    alpha = z' * w;                    
    gamma = z_old' * w;                
    w     = w - alpha * v - gamma * v_old;

    wp   = Hsolve(w);
    beta = sqrt(abs(wp' * w));

    if ~isfinite(beta) || beta <= 0
      break;
    end

    v_old = v;
    v     = w / beta;

    % next column of R in QR = T
    gammat  = conj(s_veryold) * gamma;
    deltatt = c_veryold * gamma;

    deltat  = conj(c_old) * deltatt + conj(s_old) * alpha;
    alphatt = -s_old * deltatt + c_old * alpha;

    alphat = sqrt(abs(alphatt)^2 + abs(beta)^2);

    if ~isfinite(alphatt) || alphatt == 0
      break;
    end

    c      = alphatt / alphat;
    s      = beta / alphat;

    % next column of Y and the next iterate
    y = (1 / alphatt) * (z - deltat * y_old - gammat * y_veryold);
    x = x + epsilont * y + (epsilon - epsilont_old * (alphat_old / alphatt_old)) * y_old;

    % ||r_m||_{H^-1} = beta_m |y_m(m)|, with y_m(m) = epsilont / alphatt
    est = beta * abs(epsilont / alphatt) / beta0;

    epsilont_old = epsilont;
    epsilon      = conj(c) * epsilont;
    epsilont     = -s * epsilont;

    % update for next iteration
    z_old = z;
    z     = wp / beta;

    y_veryold = y_old;
    y_old     = (alphatt / alphat) * y;

    s_veryold = s_old;  c_veryold = c_old;
    s_old     = s;      c_old     = c;

    alphat_old  = alphat;
    alphatt_old = alphatt;
  end

  iter = m - 1;
end
