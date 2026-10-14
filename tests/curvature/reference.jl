# Reference solutions for the curvature verification tests (see
# Scheme_curvature_note.tex, "Proposed numerical tests").
#
# Steady, frictionless, supercritical flow over a one-dimensional bed b(x).
# H is the bed-normal flow depth (Kestrel's Hn) and u the horizontal velocity.
# The steady equations of the family of systems considered are
#
#    (gam^k H u)_x = 0,
#    u u_x = -(g/gam^2) (b_x + (H/gam)_x) + sigma b_x b_xx u^2 / gam^2,
#
# with gam = sqrt(1 + b_x^2) and
#    (k, sigma) = (1, +1): scheme as published in Langham et al. (2023) ("scheme A"),
#    (k, sigma) = (1,  0): gamma-weighted mass and momentum, no curvature term
#                          ("geometric factors = no curvature"),
#    (k, sigma) = (0,  0): model of Langham et al. (2023),
#    (k, sigma) = (1, -1): corrected system ("geometric factors = on").
# The last also has the exact Bernoulli solution computed by valley_exact.

# Gaussian bed b = A exp(-x^2/w^2) and its first two derivatives.
struct GaussBed
   A::Float64
   w::Float64
end
bed(b::GaussBed, x) = b.A * exp(-(x / b.w)^2)
bed_x(b::GaussBed, x) = -2x / b.w^2 * bed(b, x)
bed_xx(b::GaussBed, x) = (4x^2 / b.w^4 - 2 / b.w^2) * bed(b, x)

# du/dx for the steady system (k, sigma), with H eliminated through the mass
# flux q = gam^k H u.
function steady_dudx(bd, x, u, q, k, sigma, g)
   bx = bed_x(bd, x)
   bxx = bed_xx(bd, x)
   gam = sqrt(1 + bx^2)
   gamx = bx * bxx / gam
   # H/gam = q / (gam^(k+1) u), so (H/gam)_x = -q (k+1) gamx / (gam^(k+2) u) - q u_x / (gam^(k+1) u^2)
   rhs = -(g / gam^2) * (bx - q * (k + 1) * gamx / (gam^(k + 2) * u)) + sigma * bx * bxx * u^2 / gam^2
   den = u - g * q / (gam^(k + 3) * u^2)
   den <= 0 && error("steady_dudx: flow is not supercritical at x = $x")
   return rhs / den
end

# Integrate the steady system (k, sigma) from the inflow at x0 (where H = H0,
# u = u0) and return (u, H) at the increasing points xs (all >= x0). Classical
# RK4 with steps of at most hmax.
function steady_profile(bd, xs, k, sigma; x0, H0, u0, g=9.81, hmax=1e-2)
   gam0 = sqrt(1 + bed_x(bd, x0)^2)
   q = gam0^k * H0 * u0
   f(x, u) = steady_dudx(bd, x, u, q, k, sigma, g)
   u = u0
   x = x0
   us = similar(xs, Float64)
   for (i, xt) in enumerate(xs)
      xt < x && error("steady_profile: points must be increasing and >= x0")
      n = max(1, ceil(Int, (xt - x) / hmax))
      h = (xt - x) / n
      for _ = 1:n
         k1 = f(x, u)
         k2 = f(x + h / 2, u + h / 2 * k1)
         k3 = f(x + h / 2, u + h / 2 * k2)
         k4 = f(x + h, u + h * k3)
         u += h / 6 * (k1 + 2k2 + 2k3 + k4)
         x += h
      end
      x = xt
      us[i] = u
   end
   gams = [sqrt(1 + bed_x(bd, xi)^2) for xi in xs]
   Hs = q ./ (gams .^ k .* us)
   return us, Hs
end

# Exact steady solution of the corrected system: gam H u = q and
# gam^2 u^2 / 2 + g (b + H/gam) = B, with q and B taken at the inflow x0. At each
# x the supercritical (larger) root u is found by bisection followed by Newton.
function valley_exact(bd, xs; x0, H0, u0, g=9.81)
   gam0 = sqrt(1 + bed_x(bd, x0)^2)
   q = gam0 * H0 * u0
   B = 0.5 * gam0^2 * u0^2 + g * (bed(bd, x0) + H0 / gam0)
   us = similar(xs, Float64)
   Hs = similar(xs, Float64)
   for (i, x) in enumerate(xs)
      gam = sqrt(1 + bed_x(bd, x)^2)
      b = bed(bd, x)
      F(u) = 0.5 * gam^2 * u^2 + g * (b + q / (gam^2 * u)) - B
      dF(u) = gam^2 * u - g * q / (gam^2 * u^2)
      # F is minimised at the critical speed uc; the supercritical root is above it.
      uc = cbrt(g * q / gam^4)
      F(uc) >= 0 && error("valley_exact: no supercritical solution at x = $x")
      lo, hi = uc, uc
      while F(hi) < 0
         hi *= 2
      end
      for _ = 1:60
         mid = (lo + hi) / 2
         F(mid) < 0 ? (lo = mid) : (hi = mid)
      end
      u = (lo + hi) / 2
      for _ = 1:5
         u -= F(u) / dF(u)
      end
      us[i] = u
      Hs[i] = q / (gam * u)
   end
   return us, Hs
end

# Relative L1 error of a numerical profile against a reference on the same points.
rel_l1(num, ref) = sum(abs.(num .- ref)) / sum(abs.(ref))
