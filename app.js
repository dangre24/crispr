const form = document.querySelector("#designer-form");
const sequenceInput = document.querySelector("#sequence");
const enzymeInput = document.querySelector("#enzyme");
const fileInput = document.querySelector("#sequence-file");
const mismatchInput = document.querySelector("#max-mismatches");
const mismatchValue = document.querySelector("#mismatch-value");
const loadExample = document.querySelector("#load-example");
const summary = document.querySelector("#summary");
const metrics = document.querySelector("#metrics");
const genomeViewer = document.querySelector("#genome-viewer");
const guideList = document.querySelector("#guide-list");
const exportCsv = document.querySelector("#export-csv");
let latestData = null;

const exampleSequence = [
  "TTTAACGTTGACCGTACGATCGATCGACCTGACGTTAGGCTA",
  "ATGACCGTACGTTAGCTAGCTGACCTTACGATCGGATCGAACCTGACGTTGACCGG",
  "TACGATGCTAGGCTAACCGTACCGTTACGGAATTCGATCGATCGGCTTACGTAGGCTA",
  "CCGATCGTACGATCGTACCTAGGATCGATGCTTACCGGATCGTACGTTAGCGGATCGA",
  "TTACGATCGATCGGATCGTACCGATGCTAGCTAGGCTAACGTACGATCGGATGCTAAC",
].join("");

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function metric(label, value) {
  return `<div class="metric"><strong>${escapeHtml(value)}</strong><span>${escapeHtml(label)}</span></div>`;
}

function parseSequenceFile(text) {
  if (text.trimStart().startsWith("LOCUS") || text.toUpperCase().includes("ORIGIN")) {
    return text
      .split(/ORIGIN/i)
      .pop()
      .split("//")[0]
      .replace(/[^ACGTacgt]/g, "")
      .toUpperCase();
  }
  return text
    .split(/\r?\n/)
    .filter((line) => !line.startsWith(">") && !line.startsWith(";"))
    .join("")
    .replace(/[^ACGTacgt]/g, "")
    .toUpperCase();
}

function riskClassName(label) {
  return `risk-${String(label || "low").toLowerCase()}`;
}

function mismatchHeatmap(item) {
  const positions = new Set(item.mismatchPositions || []);
  return `<div class="mismatch-map" title="Mismatch positions across the guide">${item.guide
    .split("")
    .map((base, index) => `<span class="${positions.has(index + 1) ? "mismatch" : ""}">${base}</span>`)
    .join("")}</div>`;
}

function renderOffTargets(items) {
  if (!items.length) {
    return `<p class="empty-small">No close off-targets found at this threshold.</p>`;
  }
  return `
    <ol class="off-targets">
      ${items
        .slice(0, 4)
        .map(
          (item) => `
          <li>
            ${mismatchHeatmap(item)}
            <span class="risk-badge ${riskClassName(item.riskClass)}">${escapeHtml(item.riskClass)}</span>
            ${item.mismatches} mismatch${
              item.mismatches === 1 ? "" : "es"
            } · seed ${item.seedMismatches} · ${item.risk}% risk · cut ${item.cutSite}
          </li>`
        )
        .join("")}
    </ol>`;
}

function renderGenomeViewer(data, activeGuide = data.guides[0]) {
  if (!data.guides.length) {
    genomeViewer.className = "genome-viewer empty";
    genomeViewer.innerHTML = "<p>No guides available to visualize.</p>";
    return;
  }

  const length = data.sequenceLength;
  const windowRadius = Number(document.querySelector("#zoom-range")?.value || 90);
  const center = activeGuide ? activeGuide.cutSite : Math.floor(length / 2);
  const start = Math.max(1, center - windowRadius);
  const end = Math.min(length, center + windowRadius);
  const visible = data.sequence.slice(start - 1, end);
  const activeStart = Number(activeGuide.range.split("-")[0]);
  const activeEnd = Number(activeGuide.range.split("-")[1]);
  const pamStart = Number(activeGuide.pamRange.split("-")[0]);
  const pamEnd = Number(activeGuide.pamRange.split("-")[1]);
  const deleted = activeGuide.knockout.deletedRange.split("-").map(Number);

  const sequence = visible
    .split("")
    .map((base, offset) => {
      const pos = start + offset;
      const classes = ["base"];
      if (pos >= activeStart && pos <= activeEnd) classes.push("guide-base");
      if (pos >= pamStart && pos <= pamEnd) classes.push("pam-base");
      if (pos === activeGuide.cutSite) classes.push("cut-base");
      if (pos >= deleted[0] && pos <= deleted[1]) classes.push("delete-base");
      return `<span class="${classes.join(" ")}" title="Position ${pos}">${base}</span>`;
    })
    .join("");

  const tracks = data.tracks
    .map((track) => {
      const left = Math.max(0, ((track.start - start) / Math.max(1, end - start + 1)) * 100);
      const right = Math.min(100, ((track.end - start + 1) / Math.max(1, end - start + 1)) * 100);
      const width = Math.max(2, right - left);
      return `<span class="track ${track.type}" style="left:${left}%;width:${width}%">${escapeHtml(track.label)}</span>`;
    })
    .join("");

  genomeViewer.className = "genome-viewer";
  genomeViewer.innerHTML = `
    <div class="viewer-head">
      <div>
        <h3>Interactive Genome Viewer</h3>
        <p>Viewing bases ${start}-${end}; active guide cuts at ${activeGuide.cutSite}.</p>
      </div>
      <label class="zoom-control">Zoom
        <input id="zoom-range" type="range" min="35" max="180" value="${windowRadius}" />
      </label>
    </div>
    <div class="ruler"><span>${start}</span><span>${Math.round((start + end) / 2)}</span><span>${end}</span></div>
    <div class="track-row">${tracks}</div>
    <div class="sequence-view">${sequence}</div>
    <div class="legend">
      <span><i class="guide-base"></i> guide</span>
      <span><i class="pam-base"></i> PAM</span>
      <span><i class="cut-base"></i> cut</span>
      <span><i class="delete-base"></i> deleted</span>
    </div>`;

  document.querySelector("#zoom-range").addEventListener("input", () => {
    renderGenomeViewer(data, activeGuide);
  });
}

function renderGuides(guides) {
  if (!guides.length) {
    guideList.className = "guide-list empty";
    guideList.innerHTML = "<p>No NGG PAM sites produced a full 20 nt guide.</p>";
    return;
  }

  guideList.className = "guide-list";
  guideList.innerHTML = guides
    .map(
      (guide, index) => `
      <article class="guide-card" data-guide-index="${index}">
        <div class="guide-top">
          <div>
            <div class="guide-sequence">${index + 1}. ${escapeHtml(guide.guide)}</div>
            <div class="guide-meta">
              <span class="pill">${escapeHtml(guide.enzyme)}</span>
              <span class="pill">PAM ${escapeHtml(guide.pam)}</span>
              <span class="pill">Strand ${escapeHtml(guide.strand)}</span>
              <span class="pill">Guide ${escapeHtml(guide.range)}</span>
              <span class="pill">Cut ${guide.cutSite}</span>
              <span class="pill">GC ${guide.gc}%</span>
            </div>
          </div>
          <div class="score-stack">
            <div class="score efficiency"><strong>${guide.efficiency}%</strong><span>efficiency</span></div>
            <div class="score specificity"><strong>${guide.specificity}%</strong><span>specificity</span></div>
            <div class="score success"><strong>${guide.successProbability}%</strong><span>success</span></div>
          </div>
        </div>
        <div class="guide-details">
          <div class="detail-box">
            <h3>Mutation Consequence</h3>
            <p>
              Deletes bases ${escapeHtml(guide.knockout.deletedRange)}
              (${escapeHtml(guide.knockout.deletedBases || "none")});
              length change ${guide.knockout.lengthChange} bp.
              ${escapeHtml(guide.knockout.prediction)}
            </p>
            <p class="protein-preview"><b>WT protein:</b> ${escapeHtml(guide.knockout.wildProteinPreview || "n/a")}</p>
            <p class="protein-preview"><b>Edited:</b> ${escapeHtml(guide.knockout.editedProteinPreview || "n/a")}</p>
          </div>
          <div class="detail-box">
            <h3>Off-Target Predictor</h3>
            ${renderOffTargets(guide.offTargets)}
          </div>
        </div>
      </article>`
    )
    .join("");

  document.querySelectorAll(".guide-card").forEach((card) => {
    card.addEventListener("click", () => {
      const guide = guides[Number(card.dataset.guideIndex)];
      document.querySelectorAll(".guide-card").forEach((item) => item.classList.remove("active"));
      card.classList.add("active");
      if (latestData) renderGenomeViewer(latestData, guide);
    });
  });
  document.querySelector(".guide-card")?.classList.add("active");
}

async function analyze(event) {
  event.preventDefault();
  summary.textContent = "Analyzing sequence...";
  guideList.className = "guide-list empty";
  guideList.innerHTML = "<p>Scanning PAMs and scoring guides.</p>";
  metrics.innerHTML = "";

  const response = await fetch("/api/analyze", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      sequence: sequenceInput.value,
      enzyme: enzymeInput.value,
      maxMismatches: Number(mismatchInput.value),
    }),
  });
  const data = await response.json();
  latestData = response.ok ? data : null;
  exportCsv.disabled = !response.ok || !data.guides.length;

  if (!response.ok) {
    summary.textContent = "Analysis needs a longer DNA sequence.";
    guideList.className = "guide-list";
    guideList.innerHTML = `<div class="error">${escapeHtml(data.error || "Unable to analyze sequence.")}</div>`;
    return;
  }

  summary.textContent = `Ranked ${data.guideCount} guide candidate${
    data.guideCount === 1 ? "" : "s"
  } from ${data.sequenceLength} bp.`;
  metrics.innerHTML =
    metric("Sequence length", `${data.sequenceLength} bp`) +
    metric("PAM sites", data.pamCount) +
    metric("Guides shown", data.guides.length);
  renderGenomeViewer(data);
  renderGuides(data.guides);
}

mismatchInput.addEventListener("input", () => {
  mismatchValue.textContent = mismatchInput.value;
});

loadExample.addEventListener("click", () => {
  sequenceInput.value = exampleSequence;
  sequenceInput.focus();
});

fileInput.addEventListener("change", async () => {
  const file = fileInput.files?.[0];
  if (!file) return;
  sequenceInput.value = parseSequenceFile(await file.text());
});

exportCsv.addEventListener("click", () => {
  if (!latestData?.guides?.length) return;
  const headers = [
    "guide",
    "enzyme",
    "pam",
    "strand",
    "range",
    "cutSite",
    "gc",
    "efficiency",
    "specificity",
    "successProbability",
    "offTargetCount",
    "mutationPrediction",
  ];
  const rows = latestData.guides.map((guide) =>
    headers
      .map((key) => {
        const value = key === "mutationPrediction" ? guide.knockout.prediction : guide[key];
        return `"${String(value ?? "").replaceAll('"', '""')}"`;
      })
      .join(",")
  );
  const blob = new Blob([[headers.join(","), ...rows].join("\n")], { type: "text/csv" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "crispr_guides.csv";
  link.click();
  URL.revokeObjectURL(link.href);
});

form.addEventListener("submit", analyze);
