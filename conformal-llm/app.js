/* Conformal question explorer. All computation is client-side from data.json:
   per-question option probabilities plus the calibrated thresholds q_hat
   computed offline on the calibration split. */
(function () {
  "use strict";
  var LETTERS = ["A", "B", "C", "D"];
  var DATA = null;

  function $(id) { return document.getElementById(id); }

  function scoreOf(q, letter, method) {
    if (method === "lac") return 1 - q.probs[letter];
    // APS: cumulative mass in descending probability order up to `letter`
    var order = LETTERS.slice().sort(function (a, b) { return q.probs[b] - q.probs[a]; });
    var cum = 0;
    for (var i = 0; i < order.length; i++) {
      cum += q.probs[order[i]];
      if (order[i] === letter) return cum;
    }
    return cum;
  }

  function predictionSet(q, qhat, method) {
    return LETTERS.filter(function (L) { return scoreOf(q, L, method) <= qhat; });
  }

  function shortSubject(s) { return s.replace(/_/g, " "); }

  function renderQuestion() {
    var qi = parseInt($("q-select").value, 10);
    var alpha = parseFloat($("alpha-select").value);
    var method = $("method-select").value;
    var q = DATA.questions[qi];
    var qhat = DATA.meta.q_hat[method][String(alpha)];
    var set = predictionSet(q, qhat, method);
    var inSet = set.indexOf(q.gold) !== -1;

    $("q-subject").textContent = shortSubject(q.subject) + "  ·  " + q.id;
    $("q-text").textContent = q.question;
    var opts = $("options");
    opts.innerHTML = "";
    LETTERS.forEach(function (L, i) {
      var p = q.probs[L];
      var div = document.createElement("div");
      div.className = "option" + (set.indexOf(L) !== -1 ? " in-set" : "");
      var tag = set.indexOf(L) !== -1 ? '<span class="tag">IN SET</span>' : "";
      div.innerHTML =
        '<span class="letter">' + L + (L === q.gold ? " &#10003;" : "") + "</span>" +
        '<span><span class="bar-wrap"><span class="bar" style="width:' +
        (p * 100).toFixed(1) + '%"></span></span><br>' +
        '<span class="choice-text"></span>' + tag + ' <span class="pct">' +
        (p * 100).toFixed(1) + "%</span></span>";
      div.querySelector(".choice-text").textContent = q.choices[i];
      opts.appendChild(div);
    });

    var v = $("verdict");
    v.className = inSet ? "" : "miss";
    if (set.length === 0) {
      v.textContent = "Empty set: every option scored above the calibrated threshold, " +
        "so nothing qualifies. This counts as a miss. It happens when the model is " +
        "confidently wrong, and the coverage guarantee already accounts for it.";
    } else {
      v.textContent = inSet
        ? "The correct answer (" + q.gold + ") is in the set. Set size " + set.length + " of 4."
        : "The correct answer (" + q.gold + ") is NOT in this set. " +
          "That happens about " + Math.round(alpha * 100) + "% of the time, by design.";
    }
    $("guarantee-line").textContent =
      "Threshold q\u0302 = " + qhat.toFixed(3) + " calibrated on " +
      DATA.meta.n_calibration + " held-out questions (" + method.toUpperCase() +
      "). Guarantee: a set built this way contains the truth with probability at least " +
      Math.round((1 - alpha) * 100) + "%.";
  }

  function barChart(svg, series, colors, labels) {
    // series: array of {label, values:[...]} ; grouped bars, 0..1 y
    var W = 420, H = 300, padL = 44, padB = 44, padT = 14, padR = 10;
    var NS = series.length, NB = series[0].values.length;
    var gw = (W - padL - padR) / NB, bw = Math.min(26, (gw - 14) / NS);
    var y = function (v) { return padT + (1 - Math.min(1, Math.max(0, v))) * (H - padT - padB); };
    var s = "";
    [0.2, 0.4, 0.6, 0.8, 1.0].forEach(function (t) {
      s += '<line x1="' + padL + '" y1="' + y(t) + '" x2="' + (W - padR) + '" y2="' + y(t) +
           '" stroke="#eee" stroke-dasharray="3,3"/>';
      s += '<text x="' + (padL - 6) + '" y="' + (y(t) + 4) + '" font-size="10" text-anchor="end" fill="#6b6259">' +
           Math.round(t * 100) + "%</text>";
    });
    // diagonal for coverage chart
    s += '<line x1="' + padL + '" y1="' + y(0) + '" x2="' + (W - padR) + '" y2="' + y(1.0) +
         '" stroke="#1c1a17" stroke-dasharray="5,4" stroke-width="1.5"/>';
    series.forEach(function (ser, si) {
      ser.values.forEach(function (v, bi) {
        var x = padL + bi * gw + (gw - NS * bw) / 2 + si * bw;
        s += '<rect x="' + x.toFixed(1) + '" y="' + y(v).toFixed(1) + '" width="' + bw +
             '" height="' + (y(0) - y(v)).toFixed(1) + '" fill="' + colors[si] + '"/>';
      });
    });
    labels.forEach(function (lab, bi) {
      var x = padL + bi * gw + gw / 2;
      s += '<text x="' + x.toFixed(1) + '" y="' + (H - padB + 18) + '" font-size="11" text-anchor="middle" fill="#1c1a17">' +
           lab + "</text>";
    });
    series.forEach(function (ser, si) {
      s += '<rect x="' + (W - padR - 150) + '" y="' + (padT + si * 18) + '" width="12" height="12" fill="' +
           colors[si] + '"/><text x="' + (W - padR - 134) + '" y="' + (padT + si * 18 + 10) +
           '" font-size="11" fill="#1c1a17">' + ser.name + "</text>";
    });
    svg.innerHTML = s;
  }

  function renderCharts() {
    var alphas = DATA.meta.alphas;
    var nomLabels = alphas.map(function (a) { return Math.round((1 - a) * 100) + "%"; });
    var lac = alphas.map(function (a) { return DATA.summary.methods.lac.by_alpha[String(a)].empirical_coverage; });
    var aps = alphas.map(function (a) { return DATA.summary.methods.aps.by_alpha[String(a)].empirical_coverage; });
    var naive = alphas.map(function (a) { return DATA.summary.naive_fixed_threshold_coverage[String(a)]; });
    barChart($("coverage-chart"),
      [{ name: "LAC", values: lac }, { name: "APS", values: aps }, { name: "naive", values: naive }],
      ["#2a5d8a", "#7a9e43", "#b0a89a"], nomLabels);

    var d = DATA.summary.methods.lac.by_alpha["0.1"].size_distribution;
    var d2 = DATA.summary.methods.aps.by_alpha["0.1"].size_distribution;
    var sizes = ["1", "2", "3", "4"];
    var maxN = Math.max.apply(null, sizes.map(function (k) { return Math.max(d[k] || 0, d2[k] || 0); }));
    var W = 420, H = 300, padL = 44, padB = 44, padT = 14, padR = 10;
    var gw = (W - padL - padR) / 4, bw = 30;
    var y = function (v) { return padT + (1 - v / maxN) * (H - padT - padB); };
    var s = "";
    for (var g = 0; g <= 4; g++) {
      var gv = Math.round(maxN * g / 4);
      s += '<text x="' + (padL - 6) + '" y="' + (y(gv) + 4) +
           '" font-size="10" text-anchor="end" fill="#6b6259">' + gv + "</text>";
      s += '<line x1="' + padL + '" y1="' + y(gv) + '" x2="' + (W - padR) + '" y2="' + y(gv) + '" stroke="#eee"/>';
    }
    [["LAC", d, "#2a5d8a", -1], ["APS", d2, "#7a9e43", 1]].forEach(function (cfg) {
      sizes.forEach(function (k, bi) {
        var v = cfg[1][k] || 0;
        var x = padL + bi * gw + gw / 2 + cfg[3] * (bw / 2 + 2) - bw / 2;
        s += '<rect x="' + x.toFixed(1) + '" y="' + y(v).toFixed(1) + '" width="' + bw + '" height="' +
             (y(0) - y(v)).toFixed(1) + '" fill="' + cfg[2] + '"/>';
      });
    });
    sizes.forEach(function (k, bi) {
      var x = padL + bi * gw + gw / 2;
      s += '<text x="' + x.toFixed(1) + '" y="' + (H - padB + 18) +
           '" font-size="11" text-anchor="middle">size ' + k + "</text>";
    });
    s += '<rect x="' + (W - padR - 120) + '" y="' + padT + '" width="12" height="12" fill="#2a5d8a"/>' +
         '<text x="' + (W - padR - 104) + '" y="' + (padT + 10) + '" font-size="11">LAC</text>' +
         '<rect x="' + (W - padR - 120) + '" y="' + (padT + 18) + '" width="12" height="12" fill="#7a9e43"/>' +
         '<text x="' + (W - padR - 104) + '" y="' + (padT + 28) + '" font-size="11">APS</text>';
    $("size-chart").innerHTML = s;

    var tb = document.querySelector("#subject-table tbody");
    tb.innerHTML = "";
    Object.keys(DATA.summary.methods.lac.by_alpha["0.1"].coverage_by_subject).forEach(function (subj) {
      var lc = DATA.summary.methods.lac.by_alpha["0.1"].coverage_by_subject[subj];
      var ac = DATA.summary.methods.aps.by_alpha["0.1"].coverage_by_subject[subj];
      var tr = document.createElement("tr");
      tr.innerHTML = "<td>" + shortSubject(subj) + "</td>" +
        '<td class="num">' + lc.n + "</td>" +
        '<td class="num">' + (lc.coverage * 100).toFixed(1) + "%</td>" +
        '<td class="num">' + (ac.coverage * 100).toFixed(1) + "%</td>" +
        '<td class="num">' + lc.mean_size.toFixed(2) + " / " + ac.mean_size.toFixed(2) + "</td>";
      tb.appendChild(tr);
    });
  }

  function init(meta, questions) {
    DATA = { meta: meta.meta, summary: meta.summary, questions: questions };
    $("meta-line").textContent =
      meta.meta.n_total + " MMLU questions \u00b7 " + meta.meta.model +
      " \u00b7 calibrated on " + meta.meta.n_calibration + " \u00b7 tested on " + meta.meta.n_test;
    var sel = $("q-select");
    questions.forEach(function (q, i) {
      var o = document.createElement("option");
      o.value = i;
      o.textContent = q.id + " — " + q.question.slice(0, 60) + "…";
      sel.appendChild(o);
    });
    sel.value = "0";
    ["q-select", "alpha-select", "method-select"].forEach(function (id) {
      $(id).addEventListener("change", renderQuestion);
    });
    $("random-btn").addEventListener("click", function () {
      sel.value = String(Math.floor(Math.random() * questions.length));
      renderQuestion();
    });
    renderQuestion();
    renderCharts();
  }

  var loads = [fetch("meta.json").then(function (r) { return r.json(); })];
  for (var i = 0; i < 4; i++) {
    loads.push(fetch("questions-" + i + ".json").then(function (r) { return r.json(); }));
  }
  Promise.all(loads).then(function (parts) {
    var questions = [];
    for (var i = 1; i < parts.length; i++) questions = questions.concat(parts[i]);
    init(parts[0], questions);
  }).catch(function () { $("meta-line").textContent = "could not load site data"; });
})();
