# Boxfish (Ostracioidea) Total Evidence Phylogenomics

Time-calibrated total evidence phylogeny of boxfishes (Ostraciidae + Aracanidae) with *Triodon macropterus* as outgroup, inferred using BEAST2 from molecular and morphological data with fossil calibrations.

**Final result:** `results/run4_beast_MCC_CA.tre` — 30 taxa, root age 64.8 Ma

---

## Taxa

30 total: 26 with molecular data + 4 morphology-only (*Acanthostracion guineensis*, *A. notacanthus*, *Ostracion cyanurus*, *Polyplacapros tyleri*). Outgroup: *Triodon macropterus*.

6 taxa excluded from molecular analysis (>95% missing data) but included via morphology only.

---

## Pipeline

### 1. Morphological matrix
- 120 characters, 26 taxa, digitized from Winterbottom & Tyler (1983) and Klassen (1995)
- `data/morpho/Ostracioidea_K95_WT83.nex`

### 2. Mitogenomes
- 13 complete mitogenomes downloaded from NCBI GenBank
- `data/mitogenomes/`
- Outgroup: *Triodon macropterus* (NC_009859.1)
- Reference: *Ostracion immaculatus* (NC_009865.1) — `reference/`

### 3. Mitochondrial locus extraction
- Each mitogenome aligned pairwise against reference using LASTZ
- LASTZ output MAF files: `analyses/lastz_out/`
- Gene coordinates parsed from reference GenBank annotation (`reference/*.gb`)
- 15 mito loci extracted per mitogenome: `data/loci/mito_extracted/`

### 4. Sequence data compilation
Three sources combined per locus with priority: **mitogenome extracted > longer of FToL or NCBI**

| Source | Location |
|--------|----------|
| Mitogenome-extracted loci | `data/loci/mito_extracted/` |
| Fish Tree of Life (FToL) mito loci | `data/loci/ftol_mito/` |
| Fish Tree of Life (FToL) nuclear loci | `data/loci/ftol_nuclear/` |
| NCBI GenBank mito loci (rentrez) | `data/loci/ncbi_mito/` |
| NCBI GenBank nuclear loci (rentrez) | `data/loci/ncbi_nuclear/` |
| Merged per-locus fastas | `data/loci/merged/` |

NCBI sequences >10 kb filtered (to exclude accidental full mitogenomes).

### 5. Alignment
- All 31 loci aligned with MAFFT (`--auto --adjustdirectionaccurately`) on Hoffman2 HPC
- Taxa with <50 real bases removed before alignment
- Missing taxa added back as `?` strings post-alignment (all 30 taxa present in every file)
- `_R_` prefix artefacts from MAFFT reverse-complementing removed before downstream use
- Aligned fastas: `data/loci/aligned/mito/` (15 loci) and `data/loci/aligned/nuclear/` (16 loci)

### 6. Supermatrix
- 30 taxa × 27,410 characters: 27,290 molecular (31 loci) + 120 morphological
- `analyses/beast/supermatrix/molecular_30taxa.nex`, `morpho_30taxa.nex`, `partitions.txt`

### 7. BEAST2 total evidence timetree (Run 4 — final)

**Model:**
- 31 molecular partitions: GTR+G (unlinked site models)
- 1 morphological partition: Lewis Mk (stateNumber=4, MM package)
- mtDNA partitions: UCLN relaxed clock (linked across 15 partitions)
- nuDNA partitions: UCLN relaxed clock (linked across 16 partitions)
- Morphology: strict clock
- Tree prior: Birth-Death (Gernhard 2008)

**Fossil calibrations:**

| Calibration | Fossil | Prior | Node |
|-------------|--------|-------|------|
| Root | *Ctenoplectus williamsi* | Uniform(53.7–120 Ma) | All 30 taxa |
| Total Ostraciidae | *Eolactoria sorbinii* | LogNormal offset=49.0 Ma | 29 taxa excl. *Triodon*, monophyletic |
| Crown Ostraciidae | *Oligolactoria bubiki* | LogNormal offset=30.28 Ma | 21 taxa |

**Run details:**
- 4 independent chains × 50M generations (run in 2 × 25M blocks on Hoffman2)
- logEvery=1000, storeEvery=2000
- XML: `analyses/beast/boxfish_beast_30taxa_v2.xml`
- Chains and logs not tracked in git (large files) — stored on Hoffman2 at `~/boxfish/run4/`

**Convergence:**
- 25% burnin applied in LogCombiner
- Combined ESS: posterior=385, prior=350, likelihood=907, morpho=318 — all >200 ✓
- Tracer diagnostics: `results/run4-beast-tracer.pdf`, `results/run4-beast-tracer-datatable`

**Result:**
- Root age: **64.8 Ma** (K-Pg boundary)
- Aracanidae and Ostraciidae monophyletic ✓
- *Triodon* correctly placed as outgroup (posterior=1.0) ✓
- MCC tree: `results/run4_beast_MCC_CA.tre`
- Figure: `results/beast-run-4-tree-fig.png`

---

## Repository Structure

```
Boxfish/
├── analyses/
│   ├── beast/
│   │   ├── boxfish_beast_30taxa_v2.xml   ← BEAST XML (current)
│   │   ├── prior_v2_MCC.tre              ← Prior-only run (monophyletic + uniform root)
│   │   ├── run4/                         ← Chain logs & combined trees (not on GitHub)
│   │   └── supermatrix/                  ← molecular_30taxa.nex, morpho_30taxa.nex, partitions.txt
│   ├── lastz_out/                        ← LASTZ MAF files (mitogenomes vs reference)
│   └── trees/
├── data/
│   ├── loci/
│   │   ├── aligned/                      ← ★ Final per-locus alignments (mito + nuclear)
│   │   ├── merged/                       ← Pre-alignment merged fastas
│   │   ├── mito_extracted/               ← Loci extracted from mitogenomes via LASTZ
│   │   ├── ftol_mito/ & ftol_nuclear/    ← Raw FToL source sequences
│   │   └── ncbi_mito/ & ncbi_nuclear/    ← Raw NCBI source sequences
│   ├── mitogenomes/                      ← Raw mitogenome FASTAs (13 taxa + outgroup)
│   └── morpho/
│       └── Ostracioidea_K95_WT83.nex     ← Morphological matrix
├── docs/                                 ← Species lists, coverage tables, taxonomy
├── literature/                           ← Reference PDFs
├── reference/                            ← O. immaculatus reference mitogenome (.fa + .gb)
├── results/                              ← ★ Final outputs
│   ├── run4_beast_MCC_CA.tre             ← Final MCC tree
│   ├── beast-run-4-tree-fig.png          ← Tree figure
│   ├── run4-beast-tracer.pdf             ← MCMC convergence diagnostics
│   └── run4-beast-tracer-datatable       ← ESS values per parameter
├── scripts/
│   └── boxfish_analysis_refactored.R     ← R script for tree visualization
└── archive/                              ← Superseded files from earlier runs
```

---

## Software

| Tool | Version | Use |
|------|---------|-----|
| LASTZ | 1.04 | Mitogenome pairwise alignment |
| MAFFT | 7.x | Per-locus multiple sequence alignment |
| BEAST2 | 2.7.7 | Bayesian total evidence timetree |
| LogCombiner | 2.7.7 | Chain merging + burnin |
| TreeAnnotator | 2.7.7 | MCC tree summarization |
| R (ape, ggtree) | — | Tree visualization |

Molecular analyses run on UCLA Hoffman2 HPC cluster. Job scripts available at `~/boxfish/` on Hoffman2.

---

## Key References

- Winterbottom R, Tyler JC (1983) Morphological matrix
- Klassen GJ (1995) Morphological matrix
- *Ctenoplectus williamsi* — root calibration (53.7 Ma)
- *Eolactoria sorbinii* — total Ostraciidae calibration (49.0 Ma)
- *Oligolactoria bubiki* — crown Ostraciidae calibration (30.28 Ma)
