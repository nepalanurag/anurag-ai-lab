/* Bayesian A/B testing calculator. Plain JS, no libraries.
 * Beta PDF/CDF and quadrature implemented by hand below.
 * Stats core is also require()-able under node for cross-checking. */
(function () {
  'use strict';

  /* ---- gamma and beta functions (Lanczos + continued fraction) ---- */
  var LANCZOS = [76.18009172947146, -86.50532032961677, 24.01409824083091,
                 -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];

  function gammaln(xx) {
    var y = xx, tmp = xx + 5.5;
    tmp -= (xx + 0.5) * Math.log(tmp);
    var ser = 1.000000000190015;
    for (var j = 0; j < 6; j++) { y += 1; ser += LANCZOS[j] / y; }
    return -tmp + Math.log(2.5066282746310005 * ser / xx);
  }

  function betacf(a, b, x) {
    var MAXIT = 200, EPS = 3e-14, FPMIN = 1e-300;
    var qab = a + b, qap = a + 1, qam = a - 1;
    var c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    var h = d, m, m2, aa, del;
    for (m = 1; m <= MAXIT; m++) {
      m2 = 2 * m;
      aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      del = d * c; h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }

  /* regularized incomplete beta I_x(a, b) = CDF of Beta(a, b) */
  function betaCdf(a, b, x) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    var bt = Math.exp(gammaln(a + b) - gammaln(a) - gammaln(b) +
                      a * Math.log(x) + b * Math.log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return bt * betacf(a, b, x) / a;
    return 1 - bt * betacf(b, a, 1 - x) / b;
  }

  function betaPdf(a, b, x) {
    if (x < 0 || x > 1) return 0;
    if (x === 0) { /* endpoint values are exact, not 0, when a == 1 */
      if (a === 1) return b;
      return a > 1 ? 0 : Infinity;
    }
    if (x === 1) {
      if (b === 1) return a;
      return b > 1 ? 0 : Infinity;
    }
    return Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) -
                    (gammaln(a) + gammaln(b) - gammaln(a + b)));
  }

  function simpson(f, lo, hi, n) { /* n even */
    var h = (hi - lo) / n, s = f(lo) + f(hi), i, x;
    for (i = 1; i < n; i++) {
      x = lo + i * h;
      s += f(x) * (i % 2 ? 4 : 2);
    }
    return s * h / 3;
  }

  /* P(pB > pA) = integral of f_B(y) * F_A(y) dy.
   * Panels concentrate where the B posterior has mass (mean +- 10 sd),
   * so peaked posteriors from large samples stay accurate. */
  function integrateAgainstBeta(g, a, b, nMain) {
    var m = a / (a + b);
    var s = Math.sqrt(a * b / ((a + b) * (a + b) * (a + b + 1)));
    var lo = Math.max(0, m - 10 * s), hi = Math.min(1, m + 10 * s);
    /* singular endpoints (param < 1): skip a negligible sliver so the
       grid never samples the pole; the UI already warns these are approximate */
    if (a < 1 && lo === 0) lo = 1e-12;
    if (b < 1 && hi === 1) hi = 1 - 1e-12;
    var h = function (y) { return g(y) * betaPdf(a, b, y); };
    var t = simpson(h, lo, hi, nMain);
    if (lo > 0) t += simpson(h, 0, lo, 64);
    if (hi < 1) t += simpson(h, hi, 1, 64);
    return t;
  }

  function clamp01(v) { return Math.min(1, Math.max(0, v)); }

  function pBeats(aA, bA, aB, bB) {
    return clamp01(integrateAgainstBeta(function (y) {
      return betaCdf(aA, bA, y);
    }, aB, bB, 512));
  }

  /* E[max(X - y, 0)] for X ~ Beta(aX, bX), exact via beta CDFs */
  function lossInner(y, aX, bX) {
    var mean = aX / (aX + bX);
    return mean * (1 - betaCdf(aX + 1, bX, y)) - y * (1 - betaCdf(aX, bX, y));
  }

  /* Expected loss of choosing an arm: E[max(other - chosen, 0)] */
  function expectedLoss(aA, bA, aB, bB, choose) {
    var v;
    if (choose === 'A') {
      v = integrateAgainstBeta(function (y) { return lossInner(y, aB, bB); },
                               aA, bA, 512);
    } else {
      v = integrateAgainstBeta(function (y) { return lossInner(y, aA, bA); },
                               aB, bB, 512);
    }
    return Math.max(0, v);
  }

  /* Abramowitz-Stegun erf, |err| <= 1.5e-7 */
  function erf(x) {
    var t = 1 / (1 + 0.3275911 * Math.abs(x));
    var y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t -
      0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return x < 0 ? -y : y;
  }
  function normCdf(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
  function normIsf(p) { /* Acklam's approximation */
    var a1 = -3.969683028665376e+01, a2 = 2.209460984245205e+02,
        a3 = -2.759285104469687e+02, a4 = 1.383577518672690e+02,
        a5 = -3.066479806614716e+01, a6 = 2.506628277459239e+00,
        b1 = -5.447609879822406e+01, b2 = 1.615858368580409e+02,
        b3 = -1.556989798598866e+02, b4 = 6.680131188771972e+01,
        b5 = -1.328068155288572e+01,
        c1 = -7.784894002430293e-03, c2 = -3.223964580411365e-01,
        c3 = -2.400758277161838e+00, c4 = -2.549732539343734e+00,
        c5 = 4.374664141464968e+00, c6 = 2.938163982698783e+00,
        d1 = 7.784695709041462e-03, d2 = 3.224671290700398e-01,
        d3 = 2.445134137142996e+00, d4 = 3.754408661907416e+00;
    var q, r, x;
    if (p < 0.02425) {
      q = Math.sqrt(-2 * Math.log(p));
      x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
          ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    } else if (p <= 0.97575) {
      q = p - 0.5; r = q * q;
      x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
          (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
    } else {
      q = Math.sqrt(-2 * Math.log(1 - p));
      x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
           ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    }
    return x;
  }

  /* two-sided pooled two-proportion z-test p-value */
  function freqPvalue(sA, nA, sB, nB) {
    var p1 = sA / nA, p2 = sB / nB, p = (sA + sB) / (nA + nB);
    var se = Math.sqrt(p * (1 - p) * (1 / nA + 1 / nB));
    if (se === 0) return 1;
    var z = Math.abs(p2 - p1) / se;
    return 2 * (1 - normCdf(z));
  }

  function classicalN(pA, lift) {
    var pB = pA * (1 + lift), pbar = 0.5 * (pA + pB), d = pB - pA;
    var z = normIsf(0.975) + normIsf(0.8);
    return Math.ceil(2 * pbar * (1 - pbar) * z * z / (d * d));
  }

  /* Box-Muller */
  var spare = null;
  function gauss() {
    if (spare !== null) { var v = spare; spare = null; return v; }
    var u, v2, s;
    do { u = Math.random() * 2 - 1; v2 = Math.random() * 2 - 1; s = u * u + v2 * v2; }
    while (s >= 1 || s === 0);
    var m = Math.sqrt(-2 * Math.log(s) / s);
    spare = v2 * m;
    return u * m;
  }
  /* normal approximation to Binomial(n, p); fine for planning, labeled as such */
  function simBinom(n, p) {
    var s = Math.round(n * p + Math.sqrt(n * p * (1 - p)) * gauss());
    return Math.min(n, Math.max(0, s));
  }

  var core = {
    betaPdf: betaPdf, betaCdf: betaCdf,
    pBeats: pBeats, expectedLoss: expectedLoss,
    freqPvalue: freqPvalue, classicalN: classicalN,
    simBinom: simBinom, normIsf: normIsf
  };

  /* ---- page wiring ---- */
  function val(id) { return parseFloat(document.getElementById(id).value); }
  function pct(x) { return (x * 100).toFixed(2) + '%'; }
  function el(id) { return document.getElementById(id); }

  function posterior(a0, b0, s, n) { return [a0 + s, b0 + n - s]; }

  function calculate() {
    var sA = val('sA'), nA = val('nA'), sB = val('sB'), nB = val('nB');
    var a0 = val('pa'), b0 = val('pb'), eps = val('eps');
    var err = el('err');
    err.textContent = '';
    if (!(Number.isInteger(sA) && Number.isInteger(nA) &&
          Number.isInteger(sB) && Number.isInteger(nB) &&
          sA >= 0 && nA > 0 && sB >= 0 && nB > 0 && sA <= nA && sB <= nB &&
          a0 > 0 && b0 > 0 && eps > 0)) {
      err.textContent = 'Check the inputs: successes and trials must be whole numbers, 0 <= successes <= trials, prior a,b > 0.';
      return;
    }
    var qA = posterior(a0, b0, sA, nA), qB = posterior(a0, b0, sB, nB);
    var mA = qA[0] / (qA[0] + qA[1]), mB = qB[0] / (qB[0] + qB[1]);
    var pb = pBeats(qA[0], qA[1], qB[0], qB[1]);
    var lossA = expectedLoss(qA[0], qA[1], qB[0], qB[1], 'A');
    var lossB = expectedLoss(qA[0], qA[1], qB[0], qB[1], 'B');
    var pv = freqPvalue(sA, nA, sB, nB);

    el('postA').textContent = 'Beta(' + qA[0] + ', ' + qA[1] + '), mean ' + pct(mA);
    el('postB').textContent = 'Beta(' + qB[0] + ', ' + qB[1] + '), mean ' + pct(mB);
    el('pBeatB').textContent = pct(pb);
    el('pBeatA').textContent = pct(1 - pb);
    el('lossA').textContent = pct(lossA) + ' of conversion rate';
    el('lossB').textContent = pct(lossB) + ' of conversion rate';
    el('pval').textContent = pv < 0.0001 ? pv.toExponential(2) : pv.toFixed(4);

    var dec = el('decision'), why = el('decisionWhy');
    if (lossA < eps && lossB < eps) {
      var pick = lossA <= lossB ? 'A' : 'B';
      dec.textContent = 'Either arm is fine.';
      why.textContent = 'Expected loss is under your ' + pct(eps) +
        ' threshold for both arms, so the choice barely matters. ' +
        'Pick ' + pick + ' (slightly lower loss) or the cheaper one to run.';
    } else if (lossB < eps) {
      dec.textContent = 'Ship B.';
      why.textContent = 'If B is actually worse, you lose about ' + pct(lossB) +
        ' of conversion rate on average, under your ' + pct(eps) + ' threshold.';
    } else if (lossA < eps) {
      dec.textContent = 'Keep A.';
      why.textContent = 'If A is actually worse, you lose about ' + pct(lossA) +
        ' of conversion rate on average, under your ' + pct(eps) + ' threshold.';
    } else {
      dec.textContent = 'Keep running.';
      why.textContent = 'Neither arm is safe to pick yet: expected loss is ' +
        pct(lossA) + ' for A and ' + pct(lossB) + ' for B, both above ' + pct(eps) + '.';
    }
    if (Math.min(qA[0], qA[1], qB[0], qB[1]) < 1) {
      err.textContent = 'Note: a posterior parameter is below 1, so the density ' +
        'blows up at an endpoint and the quadrature is approximate here.';
    }
    drawPosteriors(qA, qB);
  }

  function drawPosteriors(qA, qB) {
    var cv = el('plot'), ctx = cv.getContext('2d');
    var W = cv.width, H = cv.height, padL = 46, padB = 30, padT = 12, padR = 10;
    ctx.clearRect(0, 0, W, H);
    var mA = qA[0] / (qA[0] + qA[1]), mB = qB[0] / (qB[0] + qB[1]);
    var sdA = Math.sqrt(qA[0] * qA[1] / ((qA[0] + qA[1]) * (qA[0] + qA[1]) * (qA[0] + qA[1] + 1)));
    var sdB = Math.sqrt(qB[0] * qB[1] / ((qB[0] + qB[1]) * (qB[0] + qB[1]) * (qB[0] + qB[1] + 1)));
    var lo = Math.max(0, Math.min(mA - 5 * sdA, mB - 5 * sdB));
    var hi = Math.min(1, Math.max(mA + 5 * sdA, mB + 5 * sdB));
    if (hi - lo < 1e-9) { hi = lo + 1e-9; }
    var N = 240, i, x, ya, yb, maxY = 0, ysA = [], ysB = [];
    for (i = 0; i <= N; i++) {
      x = lo + (hi - lo) * i / N;
      ya = betaPdf(qA[0], qA[1], x); yb = betaPdf(qB[0], qB[1], x);
      ysA.push(ya); ysB.push(yb);
      if (ya > maxY) maxY = ya;
      if (yb > maxY) maxY = yb;
    }
    var X = function (v) { return padL + (v - lo) / (hi - lo) * (W - padL - padR); };
    var Y = function (v) { return H - padB - v / maxY * (H - padT - padB); };
    /* axes */
    ctx.strokeStyle = '#999'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(padL, H - padB); ctx.lineTo(W - padR, H - padB); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(padL, padT); ctx.lineTo(padL, H - padB); ctx.stroke();
    ctx.fillStyle = '#555'; ctx.font = '11px system-ui';
    for (i = 0; i <= 4; i++) {
      x = lo + (hi - lo) * i / 4;
      ctx.fillText((x * 100).toFixed(1) + '%', X(x) - 14, H - padB + 16);
    }
    function curve(ys, color) {
      ctx.strokeStyle = color; ctx.lineWidth = 2;
      ctx.beginPath();
      for (i = 0; i <= N; i++) {
        x = X(lo + (hi - lo) * i / N);
        if (i === 0) ctx.moveTo(x, Y(ys[i])); else ctx.lineTo(x, Y(ys[i]));
      }
      ctx.stroke();
    }
    curve(ysA, '#1f5fa8');
    curve(ysB, '#b03a2e');
    ctx.font = '12px system-ui';
    ctx.fillStyle = '#1f5fa8'; ctx.fillText('A', W - padR - 34, padT + 12);
    ctx.fillStyle = '#b03a2e'; ctx.fillText('B', W - padR - 14, padT + 12);
    ctx.fillStyle = '#555'; ctx.font = '11px system-ui';
    ctx.fillText('posterior density', padL + 4, padT + 12);
  }

  function runPlanner() {
    var pA = val('ppA'), lift = val('plift') / 100, tau = val('ptau'),
        target = val('ptarget') / 100, a0 = val('ppa'), b0 = val('ppb');
    var out = el('planOut'), err = el('planErr');
    err.textContent = '';
    if (!(pA > 0 && pA < 1 && lift > 0 && tau > 0 && tau < 1 && target > 0 && target < 1)) {
      err.textContent = 'Check the inputs: baseline in (0,1), lift > 0, thresholds in (0,1).';
      return;
    }
    out.textContent = 'running simulations...';
    setTimeout(function () {
      var pB = pA * (1 + lift);
      var nC = classicalN(pA, lift);
      var grid = [0.5, 0.75, 1, 1.25, 1.5, 2].map(function (m) {
        return Math.round(nC * m / 100) * 100;
      });
      grid = grid.filter(function (v, i) { return grid.indexOf(v) === i && v > 200; });
      var SIMS = 120, rows = [], n, s, hits, qA, qB, pb, a;
      for (var g = 0; g < grid.length; g++) {
        n = grid[g]; hits = 0;
        for (s = 0; s < SIMS; s++) {
          var xa = simBinom(n, pA), xb = simBinom(n, pB);
          qA = [a0 + xa, b0 + n - xa];
          qB = [a0 + xb, b0 + n - xb];
          pb = pBeats(qA[0], qA[1], qB[0], qB[1]);
          if (pb > tau) hits++;
        }
        a = hits / SIMS;
        rows.push({ n: n, assurance: a });
        if (a >= target) break;
      }
      var found = rows.filter(function (r) { return r.assurance >= target; })[0];
      var html = rows.map(function (r) {
        return 'n=' + r.n + '/arm: assurance ' + pct(r.assurance);
      }).join('<br>');
      if (found) {
        html += '<br><b>Plan for ' + found.n + ' visitors per arm.</b> ' +
          'About ' + pct(found.assurance) + ' of simulated trials reach ' +
          'P(B beats A) > ' + pct(tau) + ' at that size.';
      } else {
        html += '<br><b>No grid point reached ' + pct(target) +
          ' assurance.</b> Try a larger lift assumption or more traffic.';
      }
      out.innerHTML = html;
    }, 30);
  }

  if (typeof document !== 'undefined') {
    document.getElementById('calcBtn').addEventListener('click', calculate);
    document.getElementById('planBtn').addEventListener('click', runPlanner);
    calculate();
  }

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = core;
  } else if (typeof window !== 'undefined') {
    window.bayesAB = core;
  }
})();
