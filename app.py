from __future__ import annotations

import json
import math
import os
import random
import re
from dataclasses import dataclass
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import numpy as np


ROOT = Path(__file__).parent
STATIC_DIR = ROOT / "static"
BASES = "ACGT"
DNA_RE = re.compile(r"[^ACGTacgt]")


ENZYMES = {
    "spcas9": {
        "name": "SpCas9",
        "pam": "NGG",
        "guide_length": 20,
        "orientation": "downstream",
        "cut_offset": -3,
    },
    "sacas9": {
        "name": "SaCas9",
        "pam": "NNGRRT",
        "guide_length": 21,
        "orientation": "downstream",
        "cut_offset": -3,
    },
    "cas12a": {
        "name": "Cas12a",
        "pam": "TTTV",
        "guide_length": 23,
        "orientation": "upstream",
        "cut_offset": 18,
    },
}


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def clean_sequence(value: str) -> str:
    return DNA_RE.sub("", value or "").upper()


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTRYNVB", "TGCAYRNBV"))[::-1]


def iupac_match(base: str, code: str) -> bool:
    mapping = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "N": "ACGT",
        "R": "AG",
        "Y": "CT",
        "V": "ACG",
        "B": "CGT",
    }
    return base in mapping.get(code, code)


def pam_matches(seq: str, pattern: str) -> bool:
    return len(seq) == len(pattern) and all(
        iupac_match(base, code) for base, code in zip(seq, pattern)
    )


def gc_fraction(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def longest_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    best = run = 1
    for index in range(1, len(seq)):
        if seq[index] == seq[index - 1]:
            run += 1
            best = max(best, run)
        else:
            run = 1
    return best


def parse_uploaded_sequence(value: str) -> str:
    text = value or ""
    if text.lstrip().startswith("LOCUS") or "ORIGIN" in text[:1000].upper():
        origin = text.upper().split("ORIGIN", 1)[-1].split("//", 1)[0]
        return clean_sequence(origin)
    lines = []
    for line in text.splitlines():
        if line.startswith(">") or line.startswith(";"):
            continue
        lines.append(line)
    return clean_sequence("".join(lines))


def translate_dna(sequence: str) -> str:
    usable = sequence[: len(sequence) - (len(sequence) % 3)]
    return "".join(CODON_TABLE.get(usable[index : index + 3], "X") for index in range(0, len(usable), 3))


@dataclass(frozen=True)
class Guide:
    guide: str
    pam: str
    strand: str
    enzyme: str
    start: int
    end: int
    pam_start: int
    pam_end: int
    cut_site: int


class TinyLogisticModel:
    def __init__(self, weights: np.ndarray, bias: float):
        self.weights = weights
        self.bias = bias

    def predict(self, features: list[float]) -> float:
        score = float(np.dot(self.weights, np.array(features)) + self.bias)
        return 1 / (1 + math.exp(-max(min(score, 35), -35)))


def guide_features(guide: str) -> list[float]:
    gc = gc_fraction(guide)
    seed = guide[-12:] if len(guide) >= 12 else guide
    length = len(guide) or 1
    return [
        gc,
        abs(gc - 0.5),
        1.0 if 0.4 <= gc <= 0.65 else 0.0,
        guide.count("G") / length,
        guide.count("A") / length,
        seed.count("G") / len(seed) if seed else 0.0,
        1.0 if guide[-1] == "G" else 0.0,
        1.0 if "TTTT" in guide else 0.0,
        min(longest_homopolymer(guide), 6) / 6,
    ]


def off_target_features(mismatches: int, seed_mismatches: int, gc: float) -> list[float]:
    return [
        mismatches / 4,
        seed_mismatches / 12,
        1.0 if mismatches <= 2 else 0.0,
        1.0 if seed_mismatches == 0 else 0.0,
        abs(gc - 0.5),
    ]


def train_logistic(
    rows: list[list[float]], labels: list[int], epochs: int = 900, rate: float = 0.25
) -> TinyLogisticModel:
    x = np.array(rows, dtype=float)
    y = np.array(labels, dtype=float)
    weights = np.zeros(x.shape[1], dtype=float)
    bias = 0.0
    for _ in range(epochs):
        logits = np.clip(x @ weights + bias, -35, 35)
        preds = 1 / (1 + np.exp(-logits))
        error = preds - y
        weights -= rate * ((x.T @ error) / len(x) + 0.015 * weights)
        bias -= rate * float(np.mean(error))
    return TinyLogisticModel(weights, bias)


def build_efficiency_model() -> TinyLogisticModel:
    rng = random.Random(19)
    rows: list[list[float]] = []
    labels: list[int] = []
    for _ in range(900):
        guide = "".join(rng.choice(BASES) for _ in range(20))
        gc = gc_fraction(guide)
        score = 0.0
        score += 2.0 if 0.42 <= gc <= 0.62 else -1.0
        score += 0.7 if guide[19] == "G" else 0.0
        score += 0.35 if guide[-12:].count("G") >= 3 else 0.0
        score -= 1.5 if "TTTT" in guide else 0.0
        score -= 1.1 if longest_homopolymer(guide) >= 5 else 0.0
        score -= abs(gc - 0.52) * 2.6
        score += rng.uniform(-0.55, 0.55)
        rows.append(guide_features(guide))
        labels.append(1 if score > 0.35 else 0)
    return train_logistic(rows, labels)


def build_off_target_model() -> TinyLogisticModel:
    rows: list[list[float]] = []
    labels: list[int] = []
    for mismatches in range(0, 5):
        for seed_mismatches in range(0, min(12, mismatches) + 1):
            for gc in [0.3, 0.4, 0.5, 0.6, 0.7]:
                risk = 0.0
                risk += 2.5 if mismatches <= 1 else 0.0
                risk += 1.2 if mismatches == 2 else 0.0
                risk += 1.4 if seed_mismatches == 0 else 0.0
                risk -= 1.0 if seed_mismatches >= 3 else 0.0
                risk -= mismatches * 0.75
                risk -= abs(gc - 0.5) * 0.8
                rows.append(off_target_features(mismatches, seed_mismatches, gc))
                labels.append(1 if risk > 0.6 else 0)
    return train_logistic(rows, labels, epochs=700, rate=0.18)


EFFICIENCY_MODEL = build_efficiency_model()
OFF_TARGET_MODEL = build_off_target_model()


def find_guides(sequence: str, enzyme_key: str = "spcas9") -> list[Guide]:
    enzyme = ENZYMES.get(enzyme_key, ENZYMES["spcas9"])
    guide_length = int(enzyme["guide_length"])
    pam_pattern = str(enzyme["pam"])
    pam_length = len(pam_pattern)
    guides: list[Guide] = []

    if enzyme["orientation"] == "downstream":
        for index in range(guide_length, len(sequence) - pam_length + 1):
            pam = sequence[index : index + pam_length]
            if not pam_matches(pam, pam_pattern):
                continue
            guides.append(
                Guide(
                    guide=sequence[index - guide_length : index],
                    pam=pam,
                    strand="+",
                    enzyme=str(enzyme["name"]),
                    start=index - guide_length + 1,
                    end=index,
                    pam_start=index + 1,
                    pam_end=index + pam_length,
                    cut_site=index + int(enzyme["cut_offset"]),
                )
            )

        reverse_pam = reverse_complement(pam_pattern)
        for index in range(0, len(sequence) - pam_length - guide_length + 1):
            pam_genomic = sequence[index : index + pam_length]
            if not pam_matches(pam_genomic, reverse_pam):
                continue
            protospacer = sequence[index + pam_length : index + pam_length + guide_length]
            guides.append(
                Guide(
                    guide=reverse_complement(protospacer),
                    pam=reverse_complement(pam_genomic),
                    strand="-",
                    enzyme=str(enzyme["name"]),
                    start=index + pam_length + 1,
                    end=index + pam_length + guide_length,
                    pam_start=index + 1,
                    pam_end=index + pam_length,
                    cut_site=index + pam_length - int(enzyme["cut_offset"]),
                )
            )

    if enzyme["orientation"] == "upstream":
        for index in range(0, len(sequence) - pam_length - guide_length + 1):
            pam = sequence[index : index + pam_length]
            if not pam_matches(pam, pam_pattern):
                continue
            guide_seq = sequence[index + pam_length : index + pam_length + guide_length]
            guides.append(
                Guide(
                    guide=guide_seq,
                    pam=pam,
                    strand="+",
                    enzyme=str(enzyme["name"]),
                    start=index + pam_length + 1,
                    end=index + pam_length + guide_length,
                    pam_start=index + 1,
                    pam_end=index + pam_length,
                    cut_site=index + pam_length + int(enzyme["cut_offset"]),
                )
            )

        reverse_pam = reverse_complement(pam_pattern)
        for index in range(guide_length, len(sequence) - pam_length + 1):
            pam_genomic = sequence[index : index + pam_length]
            if not pam_matches(pam_genomic, reverse_pam):
                continue
            protospacer = sequence[index - guide_length : index]
            guides.append(
                Guide(
                    guide=reverse_complement(protospacer),
                    pam=reverse_complement(pam_genomic),
                    strand="-",
                    enzyme=str(enzyme["name"]),
                    start=index - guide_length + 1,
                    end=index,
                    pam_start=index + 1,
                    pam_end=index + pam_length,
                    cut_site=index + pam_length - int(enzyme["cut_offset"]),
                )
            )
    return sorted(guides, key=lambda g: (g.cut_site, g.strand))


def mismatch_report(a: str, b: str) -> tuple[int, int, list[int]]:
    positions = [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]
    seed_mismatches = sum(1 for pos in positions if pos >= 9)
    return len(positions), seed_mismatches, positions


def risk_class(risk: float) -> str:
    if risk >= 72:
        return "High"
    if risk >= 45:
        return "Moderate"
    return "Low"


def scan_off_targets(
    sequence: str, selected: Guide, max_mismatches: int, enzyme_key: str
) -> list[dict[str, Any]]:
    sites: list[dict[str, Any]] = []
    candidates = find_guides(sequence, enzyme_key)
    for candidate in candidates:
        if (
            candidate.start == selected.start
            and candidate.end == selected.end
            and candidate.strand == selected.strand
        ):
            continue
        mismatches, seed_mismatches, positions = mismatch_report(selected.guide, candidate.guide)
        if mismatches <= max_mismatches:
            risk = OFF_TARGET_MODEL.predict(
                off_target_features(mismatches, seed_mismatches, gc_fraction(candidate.guide))
            )
            risk_percent = round(risk * 100, 1)
            sites.append(
                {
                    "guide": candidate.guide,
                    "pam": candidate.pam,
                    "strand": candidate.strand,
                    "range": f"{candidate.start}-{candidate.end}",
                    "cutSite": candidate.cut_site,
                    "mismatches": mismatches,
                    "seedMismatches": seed_mismatches,
                    "mismatchPositions": positions,
                    "risk": risk_percent,
                    "riskClass": risk_class(risk_percent),
                }
            )
    return sorted(sites, key=lambda item: (item["mismatches"], -item["risk"], item["cutSite"]))[:12]


def simulate_knockout(sequence: str, cut_site: int, deletion_size: int = 8) -> dict[str, Any]:
    left = max(0, cut_site - deletion_size // 2)
    right = min(len(sequence), left + deletion_size)
    edited = sequence[:left] + sequence[right:]
    wild_protein = translate_dna(sequence)
    edited_protein = translate_dna(edited)
    frameshift = abs(len(edited) - len(sequence)) % 3 != 0
    stop_position = edited_protein.find("*") + 1 if "*" in edited_protein else None
    return {
        "cutSite": cut_site,
        "deletedRange": f"{left + 1}-{right}",
        "deletedBases": sequence[left:right],
        "editedPreview": edited[max(0, left - 35) : min(len(edited), left + 35)],
        "lengthChange": len(edited) - len(sequence),
        "frameshift": frameshift,
        "prediction": "Frameshift mutation -> likely knockout"
        if frameshift
        else "In-frame deletion -> protein may remain partially functional",
        "wildProteinPreview": wild_protein[:80],
        "editedProteinPreview": edited_protein[:80],
        "prematureStop": stop_position,
    }


def build_exon_tracks(length: int) -> list[dict[str, Any]]:
    if length < 60:
        return [{"type": "exon", "start": 1, "end": length, "label": "exon 1"}]
    exon_one_end = max(24, int(length * 0.28))
    exon_two_start = max(exon_one_end + 12, int(length * 0.48))
    exon_two_end = min(length, int(length * 0.72))
    exon_three_start = min(length, max(exon_two_end + 12, int(length * 0.84)))
    tracks = [
        {"type": "exon", "start": 1, "end": exon_one_end, "label": "exon 1"},
        {"type": "intron", "start": exon_one_end + 1, "end": exon_two_start - 1, "label": "intron"},
        {"type": "exon", "start": exon_two_start, "end": exon_two_end, "label": "exon 2"},
    ]
    if exon_three_start < length:
        tracks.append(
            {"type": "intron", "start": exon_two_end + 1, "end": exon_three_start - 1, "label": "intron"}
        )
        tracks.append({"type": "exon", "start": exon_three_start, "end": length, "label": "exon 3"})
    return tracks


def analyze(payload: dict[str, Any]) -> dict[str, Any]:
    sequence = parse_uploaded_sequence(payload.get("sequence", ""))
    enzyme_key = str(payload.get("enzyme", "spcas9")).lower()
    enzyme = ENZYMES.get(enzyme_key, ENZYMES["spcas9"])
    max_mismatches = int(payload.get("maxMismatches", 3))
    max_mismatches = max(0, min(max_mismatches, 4))
    min_length = int(enzyme["guide_length"]) + len(str(enzyme["pam"]))
    if len(sequence) < min_length:
        return {
            "error": f"Enter at least {min_length} DNA bases so the tool can find a guide plus PAM."
        }

    guides = find_guides(sequence, enzyme_key)
    guide_rows = []
    for guide in guides:
        efficiency = EFFICIENCY_MODEL.predict(guide_features(guide.guide))
        off_targets = scan_off_targets(sequence, guide, max_mismatches, enzyme_key)
        specificity = max(0.0, 100 - sum(item["risk"] for item in off_targets[:5]) * 0.22)
        success_probability = max(
            0.0, min(100.0, (efficiency * 100 * 0.62) + (specificity * 0.32) + 6)
        )
        guide_rows.append(
            {
                "guide": guide.guide,
                "pam": guide.pam,
                "strand": guide.strand,
                "enzyme": guide.enzyme,
                "range": f"{guide.start}-{guide.end}",
                "pamRange": f"{guide.pam_start}-{guide.pam_end}",
                "cutSite": guide.cut_site,
                "gc": round(gc_fraction(guide.guide) * 100, 1),
                "efficiency": round(efficiency * 100, 1),
                "specificity": round(specificity, 1),
                "successProbability": round(success_probability, 1),
                "offTargetCount": len(off_targets),
                "offTargets": off_targets,
                "knockout": simulate_knockout(sequence, guide.cut_site),
            }
        )

    guide_rows.sort(key=lambda item: (item["efficiency"] + item["specificity"]) / 2, reverse=True)
    pam_count = sum(1 for guide in guides if guide.strand == "+") + sum(
        1 for guide in guides if guide.strand == "-"
    )
    return {
        "sequenceLength": len(sequence),
        "sequence": sequence,
        "enzyme": enzyme,
        "pamCount": pam_count,
        "guideCount": len(guide_rows),
        "guides": guide_rows[:30],
        "tracks": build_exon_tracks(len(sequence)),
        "model": {
            "efficiencyFeatures": [
                "GC balance",
                "seed composition",
                "terminal G",
                "poly-T penalty",
                "homopolymer penalty",
                "guide length normalization",
            ],
            "offTargetFeatures": ["total mismatches", "seed mismatches", "GC balance"],
        },
    }


class Handler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:
    parsed = urlparse(self.path)
    if parsed.path in {"/", "/index.html"}:
        self.send_file(ROOT / "index.html", "text/html; charset=utf-8")
    elif parsed.path in {"/styles.css", "/static/styles.css"}:
        self.send_file(ROOT / "styles.css", "text/css; charset=utf-8")
    elif parsed.path in {"/app.js", "/static/app.js"}:
        self.send_file(ROOT / "app.js", "application/javascript; charset=utf-8")
    else:
        self.send_error(HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:
        if urlparse(self.path).path != "/api/analyze":
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        length = int(self.headers.get("Content-Length", "0"))
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
            result = analyze(payload)
            self.send_json(result, HTTPStatus.BAD_REQUEST if "error" in result else HTTPStatus.OK)
        except Exception as exc:
            self.send_json({"error": str(exc)}, HTTPStatus.INTERNAL_SERVER_ERROR)

    def send_file(self, path: Path, content_type: str) -> None:
        if not path.exists():
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        body = path.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_json(self, payload: dict[str, Any], status: HTTPStatus) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, format: str, *args: Any) -> None:
        return


def main() -> None:
    port = int(os.environ.get("PORT", "8000"))
    server = ThreadingHTTPServer(("0.0.0.0", port), Handler)
    print(f"CRISPR Guide Designer running at http://127.0.0.1:{port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
