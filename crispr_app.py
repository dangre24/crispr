from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from http import HTTPStatus
import json, re, math, random

DNA_RE = re.compile(r"[^ACGTacgt]")
BASES = "ACGT"

HTML = r"""
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>CRISPR Guide RNA Designer</title>
  <style>
    body { margin:0; font-family:Arial, sans-serif; background:#f4f7f2; color:#14211f; }
    main { max-width:1150px; margin:auto; padding:32px 18px; }
    h1 { font-size:48px; margin:0 0 8px; }
    .lede { color:#5d6a67; max-width:760px; line-height:1.5; }
    .grid { display:grid; grid-template-columns:0.8fr 1.2fr; gap:18px; align-items:start; }
    .panel { background:white; border:1px solid #d8ded8; border-radius:8px; padding:18px; box-shadow:0 14px 35px #0001; }
    textarea { width:100%; min-height:300px; box-sizing:border-box; border:1px solid #c8d1cc; border-radius:8px; padding:12px; font-family:monospace; }
    button { background:#0f766e; color:white; border:0; border-radius:8px; padding:12px 16px; font-weight:700; cursor:pointer; }
    button.secondary { background:white; color:#0f766e; border:1px solid #d8ded8; }
    .controls { display:flex; gap:14px; align-items:center; margin-top:14px; flex-wrap:wrap; }
    .metrics { display:grid; grid-template-columns:repeat(3,1fr); gap:10px; margin:14px 0; }
    .metric, .card { border:1px solid #d8ded8; border-radius:8px; background:#fbfcfa; padding:12px; }
    .metric strong { display:block; font-size:24px; }
    .card { background:white; margin-bottom:12px; border-left:5px solid #0f766e; }
    .seq { font-family:monospace; font-weight:800; overflow-wrap:anywhere; }
    .pill { display:inline-block; background:#eef4f1; padding:5px 9px; border-radius:999px; margin:6px 4px 0 0; font-size:13px; font-weight:700; }
    .scores { display:flex; gap:8px; margin-top:10px; }
    .score { color:white; border-radius:8px; padding:8px 10px; text-align:center; font-weight:800; }
    .eff { background:#0f766e; } .spec { background:#2563eb; }
    .small { color:#5d6a67; font-size:14px; line-height:1.45; }
    @media(max-width:850px){ .grid,.metrics{grid-template-columns:1fr;} h1{font-size:34px;} }
  </style>
</head>
<body>
<main>
  <h1>CRISPR Guide RNA Designer</h1>
  <p class="lede">Paste a DNA sequence to find SpCas9 NGG PAM sites, generate guide RNAs, score editing efficiency, predict off-target risk, and preview a small knockout deletion.</p>

  <section class="grid">
    <div class="panel">
      <h2>Input Sequence</h2>
      <p class="small">Only A, C, G, and T are used. Spaces and numbers are ignored.</p>
      <textarea id="sequence" placeholder="Paste DNA sequence here..."></textarea>
      <div class="controls">
        <label>Max mismatches: <b id="mmText">3</b>
          <input id="mm" type="range" min="0" max="4" value="3">
        </label>
        <button type="button" onclick="loadExample()" class="secondary">Load example</button>
        <button type="button" onclick="analyze()">Analyze guides</button>
      </div>
      <p class="small">Educational prototype only. Not clinical or production genome-editing guidance.</p>
    </div>

    <div class="panel">
      <h2>Results</h2>
      <p id="summary" class="small">Run an analysis to see guide candidates.</p>
      <div id="metrics" class="metrics"></div>
      <div id="results"></div>
    </div>
  </section>
</main>

<script>
const example = [
"ATGACCGTACGTTAGCTAGCTGACCTTACGATCGGATCGAACCTGACGTTGACCGG",
"TACGATGCTAGGCTAACCGTACCGTTACGGAATTCGATCGATCGGCTTACGTAGGCTA",
"CCGATCGTACGATCGTACCTAGGATCGATGCTTACCGGATCGTACGTTAGCGGATCGA",
"TTACGATCGATCGGATCGTACCGATGCTAGCTAGGCTAACGTACGATCGGATGCTAAC"
].join("");

document.getElementById("mm").oninput = e => document.getElementById("mmText").textContent = e.target.value;

function loadExample() {
  document.getElementById("sequence").value = example;
}

function esc(x) {
  return String(x).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll(">","&gt;");
}

async function analyze() {
  const res = await fetch("/api/analyze", {
    method: "POST",
    headers: {"Content-Type":"application/json"},
    body: JSON.stringify({
      sequence: document.getElementById("sequence").value,
      maxMismatches: Number(document.getElementById("mm").value)
    })
  });
  const data = await res.json();

  if (!res.ok) {
    document.getElementById("summary").textContent = data.error || "Could not analyze sequence.";
    return;
  }

  document.getElementById("summary").textContent =
    `Ranked ${data.guideCount} guide candidates from ${data.sequenceLength} bp.`;

  document.getElementById("metrics").innerHTML = `
    <div class="metric"><strong>${data.sequenceLength} bp</strong><span>Sequence length</span></div>
    <div class="metric"><strong>${data.pamCount}</strong><span>PAM sites</span></div>
    <div class="metric"><strong>${data.guides.length}</strong><span>Guides shown</span></div>
  `;

  document.getElementById("results").innerHTML = data.guides.map((g, i) => `
    <div class="card">
      <div class="seq">${i + 1}. ${esc(g.guide)}</div>
      <span class="pill">PAM ${esc(g.pam)}</span>
      <span class="pill">Strand ${esc(g.strand)}</span>
      <span class="pill">Guide ${esc(g.range)}</span>
      <span class="pill">Cut ${g.cutSite}</span>
      <span class="pill">GC ${g.gc}%</span>
      <div class="scores">
        <div class="score eff">${g.efficiency}%<br><small>efficiency</small></div>
        <div class="score spec">${g.specificity}%<br><small>specificity</small></div>
      </div>
      <p class="small"><b>Knockout:</b> deletes bases ${esc(g.knockout.deletedRange)}
      (${esc(g.knockout.deletedBases)}); length change ${g.knockout.lengthChange} bp.</p>
      <p class="small"><b>Off-targets:</b> ${
        g.offTargets.length
          ? g.offTargets.slice(0,4).map(o => `${esc(o.guide)} · ${o.mismatches} mismatches · ${o.risk}% risk`).join("<br>")
          : "No close off-targets found at this threshold."
      }</p>
    </div>
  `).join("");
}
</script>
</body>
</html>
"""

def clean_sequence(seq):
    return DNA_RE.sub("", seq or "").upper()

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def gc_fraction(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0

def longest_homopolymer(seq):
    best = run = 1 if seq else 0
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            best = max(best, run)
        else:
            run = 1
    return best

def sigmoid(x):
    return 1 / (1 + math.exp(-max(min(x, 35), -35)))

def efficiency_score(guide):
    gc = gc_fraction(guide)
    score = 0
    score += 2.0 if 0.42 <= gc <= 0.62 else -1.0
    score += 0.7 if guide[-1] == "G" else 0
    score += 0.35 if guide[-12:].count("G") >= 3 else 0
    score -= 1.5 if "TTTT" in guide else 0
    score -= 1.1 if longest_homopolymer(guide) >= 5 else 0
    score -= abs(gc - 0.52) * 2.6
    return round(sigmoid(score) * 100, 1)

def off_target_risk(mismatches, seed_mismatches, gc):
    score = 0
    score += 2.5 if mismatches <= 1 else 0
    score += 1.2 if mismatches == 2 else 0
    score += 1.4 if seed_mismatches == 0 else 0
    score -= 1.0 if seed_mismatches >= 3 else 0
    score -= mismatches * 0.75
    score -= abs(gc - 0.5) * 0.8
    return round(sigmoid(score) * 100, 1)

def find_guides(sequence):
    guides = []

    # Plus strand: 20 nt guide followed by NGG PAM.
    for i in range(20, len(sequence) - 2):
        pam = sequence[i:i+3]
        if pam[1:] == "GG":
            guides.append({
                "guide": sequence[i-20:i],
                "pam": pam,
                "strand": "+",
                "start": i - 19,
                "end": i,
                "pam_start": i + 1,
                "pam_end": i + 3,
                "cut": i - 3
            })

    # Minus strand: genomic CCN is reverse-complement NGG.
    for i in range(0, len(sequence) - 22):
        pam_genomic = sequence[i:i+3]
        if pam_genomic[:2] == "CC":
            protospacer = sequence[i+3:i+23]
            guides.append({
                "guide": reverse_complement(protospacer),
                "pam": reverse_complement(pam_genomic),
                "strand": "-",
                "start": i + 4,
                "end": i + 23,
                "pam_start": i + 1,
                "pam_end": i + 3,
                "cut": i + 6
            })

    return sorted(guides, key=lambda g: (g["cut"], g["strand"]))

def mismatch_report(a, b):
    positions = [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]
    seed_mismatches = sum(1 for p in positions if p >= 9)
    return len(positions), seed_mismatches, positions

def simulate_knockout(sequence, cut_site, deletion_size=8):
    left = max(0, cut_site - deletion_size // 2)
    right = min(len(sequence), left + deletion_size)
    edited = sequence[:left] + sequence[right:]
    return {
        "cutSite": cut_site,
        "deletedRange": f"{left + 1}-{right}",
        "deletedBases": sequence[left:right],
        "editedPreview": edited[max(0, left - 35):min(len(edited), left + 35)],
        "lengthChange": len(edited) - len(sequence)
    }

def scan_off_targets(sequence, selected, all_guides, max_mismatches):
    hits = []
    for cand in all_guides:
        if cand["start"] == selected["start"] and cand["end"] == selected["end"] and cand["strand"] == selected["strand"]:
            continue
        mismatches, seed_mismatches, positions = mismatch_report(selected["guide"], cand["guide"])
        if mismatches <= max_mismatches:
            hits.append({
                "guide": cand["guide"],
                "pam": cand["pam"],
                "strand": cand["strand"],
                "range": f'{cand["start"]}-{cand["end"]}',
                "cutSite": cand["cut"],
                "mismatches": mismatches,
                "seedMismatches": seed_mismatches,
                "mismatchPositions": positions,
                "risk": off_target_risk(mismatches, seed_mismatches, gc_fraction(cand["guide"]))
            })
    return sorted(hits, key=lambda x: (x["mismatches"], -x["risk"], x["cutSite"]))[:12]

def analyze(payload):
    sequence = clean_sequence(payload.get("sequence", ""))
    max_mismatches = max(0, min(int(payload.get("maxMismatches", 3)), 4))

    if len(sequence) < 23:
        return {"error": "Enter at least 23 DNA bases so the tool can find a 20 nt guide plus PAM."}, HTTPStatus.BAD_REQUEST

    guides = find_guides(sequence)
    rows = []

    for guide in guides:
        off_targets = scan_off_targets(sequence, guide, guides, max_mismatches)
        specificity = max(0, 100 - sum(o["risk"] for o in off_targets[:5]) * 0.22)

        rows.append({
            "guide": guide["guide"],
            "pam": guide["pam"],
            "strand": guide["strand"],
            "range": f'{guide["start"]}-{guide["end"]}',
            "pamRange": f'{guide["pam_start"]}-{guide["pam_end"]}',
            "cutSite": guide["cut"],
            "gc": round(gc_fraction(guide["guide"]) * 100, 1),
            "efficiency": efficiency_score(guide["guide"]),
            "specificity": round(specificity, 1),
            "offTargetCount": len(off_targets),
            "offTargets": off_targets,
            "knockout": simulate_knockout(sequence, guide["cut"])
        })

    rows.sort(key=lambda g: (g["efficiency"] + g["specificity"]) / 2, reverse=True)

    return {
        "sequenceLength": len(sequence),
        "pamCount": len(guides),
        "guideCount": len(rows),
        "guides": rows[:30]
    }, HTTPStatus.OK


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path in ["/", "/index.html"]:
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(HTML.encode("utf-8"))
        else:
            self.send_error(HTTPStatus.NOT_FOUND)

    def do_POST(self):
        if self.path != "/api/analyze":
            self.send_error(HTTPStatus.NOT_FOUND)
            return

        length = int(self.headers.get("Content-Length", "0"))
        payload = json.loads(self.rfile.read(length) or b"{}")
        result, status = analyze(payload)

        body = json.dumps(result).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


import os

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8000))
    print(f"CRISPR Guide RNA Designer running on port {port}")
    ThreadingHTTPServer(("0.0.0.0", port), Handler).serve_forever()
   
