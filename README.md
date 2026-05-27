# CRISPR Guide RNA Designer + Off-Target Predictor

A lightweight Python web tool for exploring CRISPR-Cas9 guide design.

## What it does

- Accepts a DNA sequence.
- Accepts raw DNA, FASTA, GenBank-style text, and `.txt` uploads.
- Detects PAM sites on both strands for SpCas9, SaCas9, and Cas12a.
- Generates 20 nt guide RNAs.
- Shows an interactive mini genome viewer with PAM, guide, cut, deletion, and exon/intron tracks.
- Predicts cut sites and simulates a small knockout deletion.
- Scores guide efficiency with a tiny NumPy-trained logistic model.
- Predicts off-target risk by scanning near-matching guides in the submitted sequence.
- Shows mismatch heatmaps, protein translation previews, frameshift status, and CSV export.

The model is intentionally small and educational. It is useful for demonstrating biotech software ideas, not for clinical or production genome-editing decisions.

## Run it

```bash
python3 app.py
```

Then open:

```text
http://127.0.0.1:8000
```

## Files

- `app.py`: Python server, CRISPR logic, and lightweight ML scoring.
- `static/index.html`: App layout.
- `static/styles.css`: Visual design.
- `static/app.js`: Browser-side interaction.
- `requirements.txt`: Python package list for deployment.
