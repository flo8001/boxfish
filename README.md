What have we done so far?
1. Inspected morphological matrix digitized by Santini, pulled from Winterbottom & Tyler 1983 and Klassen 1995 `/morpho/Ostracioidea_K95_WT83.nex`
2. Manually downloaded mitogenomes from NCBI Genbank `/mitogenomes`
	1. aligned to reference mitogenome = `Ostracion_immaculatus_Ref_mitogenome_NC_009865.1` using LASTZ `/mitogenomes/all_mitogenomes_aligned.fasta`
	2. outgroup = `Triodon_macropterus_mitogenome_NC_009859.1`
3. Downloaded FToL mtDNA and nuDNA loci data `/ftol`
	1. extracted loci using partition coordinates `/loci/ftol_mito`, `/loci/ftol_nuclear`
4. Conducted rentrez search on NCBI for all mtDNA and nuDNA loci data
	1. stored in `/loci/ncbi_mito`, `/loci/ncbi_nuclear`
5. Extracted individual mitochondrial loci from each mitogenome using LASTZ pairwise alignments against _O. immaculatus_ reference
    1. Gene coordinates parsed from reference GenBank file `/reference/Ostracion_immaculatus_Ref_mitogenome_NC_009865.1.gb`
    2. Extracted sequences stored in `/loci/mito_extracted/`
6. Merged all three data sources (mitogenome extracted, FToL, NCBI) into per-locus fastas with priority: mitogenome > longest of FToL or NCBI
    1. Rich headers added to each sequence recording source, accession number, and sequence length
    2. Merged fastas stored in `/loci/merged/mito/` and `/loci/merged/nuclear/`
7. Aligned merged per-locus fastas using MAFFT (`--auto --adjustdirectionaccurately`) on Hoffman2
    1. Taxa with <50 real bases filtered out before alignment
    2. Missing taxa added back as `?` strings after alignment to ensure all 26 taxa present in every locus file
    3. Aligned fastas stored in `/loci/aligned/mito/` and `/loci/aligned/nuclear/`
8. Built total evidence supermatrix concatenating all 27 aligned loci + morphological matrix
    1. 26 taxa × 27,410 characters (27,290 molecular + 120 morphological)
    2. Mean missing data 63.5%
    3. Stored in `/supermatrix/`
9. Per-locus maximum likelihood gene trees inferred with IQ-TREE2 (GTR+G, 1000 ultrafast bootstraps) on Hoffman2
    1. 27/31 loci produced trees (4c4, ficd, hoxc6a, kiaa1239 skipped — insufficient taxa)
    2. Gene trees stored in `/gene_trees/`
10. Coalescent species tree inferred with ASTRAL from 27 gene trees
    1. `_R_` prefix artefacts from MAFFT reverse-complementing removed before running
    2. Stored in `/trees/astral_species_tree_final.tre`
11. Total evidence Bayesian timetree inferred with BEAST2 on Hoffman2
    1. 31 molecular partitions (GTR+G) + morphological partition (Lewis Mk)
    2. Separate UCLN relaxed clocks for mtDNA and nuDNA, strict clock for morphology
    3. Birth Death tree prior
    4. 3 fossil calibrations: _Ctenoplectus williamsi_ (53.7 Ma, root), _Eolactoria sorbinii_ (49 Ma, total Ostraciidae), _Oligolactoria bubiki_ (30.28 Ma, crown Ostraciidae)
    5. 4 independent chains × 20M generations, 25% burnin
    6. ESS 66-123 — below 200 threshold, longer runs needed
    7. MCC tree summarised with TreeAnnotator, stored in `/beast_MCC.tre`

Boxfish/
│
├── archive/                          ← Old/superseded work, kept for reference
│   ├── _FScript_TINA.R               ← Early R script (draft version)
│   ├── FScript_TINA.R                ← Early R script (cleaner version)
│   ├── modifying branch lengths.R    ← R script for manual tree editing
│   ├── comparative_boxfish_phylogeny.phy  ← Old phylip alignment
│   ├── timeCalibratedFull.phy        ← Old time-calibrated tree
│   ├── USE_THIS_TINA_timeCalibratedFull.phy  ← Final version of above
│   ├── Ostracioidea_K95_WT83.nex     ← Old morphological nexus matrix
│   ├── Ostracioidea_K95_WT83.nex copy
│   ├── Ostracioidea_K95_WT83.nex.tree
│   ├── muscle_results/               ← Old MUSCLE alignment run outputs
│   │   └── muscle-I20260206-*.*      ← Alignment, pairwise identity, tree files
│   └── individual_seqs/              ← Manually downloaded individual NCBI sequences
│       └── Genus_species_LOCUS_ACCESSION.fasta
│
├── beast/                            ← BEAST2 Bayesian analysis files
│   ├── boxfish_beast.xml             ← BEAST input: model, partitions, priors, taxa
│   ├── testbeast.xml                 ← Test/trial XML before final run
│   ├── chain_1-4_$(boxfish).log      ← MCMC parameter logs per chain (~140MB each, not on GitHub)
│   ├── chain_1-4_$(boxfish).trees    ← Posterior tree samples per chain (not on GitHub)
│   ├── combined.trees                ← All 4 chains merged (not on GitHub)
│   ├── combined_25.trees             ← 25% burnin removed, ready for analysis (not on GitHub)
│   ├── beast_MCC.tre                 ← ★ Final MCC consensus tree (main result)
│   └── boxfish_total.con.tre         ← Total evidence consensus tree
│
├── docs/                             ← Project documents and metadata
│   ├── boxfish_sp_list.xlsx          ← Master species list
│   ├── ~$boxfish_sp_list.xlsx        ← Excel temp file (can delete)
│   ├── overview_species_drag_dataset.xlsx  ← Data coverage overview per species
│   ├── ncbi_coverage.csv             ← NCBI sequence availability per species/locus
│   ├── taxonomy_lookup.csv           ← Taxonomy reference table
│   └── Boxfish Phylogenetic Methods paper copy.pdf  ← Reference paper
│
├── ftol/                             ← Fish Tree of Life sequences & IQ-TREE analysis
│   ├── Ostraciidae.nex               ← Boxfish subset from FTOL (nexus)
│   ├── Ostraciidae.phylip            ← Same, phylip format
│   ├── Ostraciidae_with_triodon.phylip  ← With outgroup added
│   ├── final_alignment.phylip        ← Final FTOL alignment used (~268MB, not on GitHub)
│   ├── final_alignment.partitions    ← Partition scheme for FTOL alignment
│   └── nuclear_iqtree.treefile       ← ★ IQ-TREE result from nuclear-only FTOL data
│
├── lastz_out/                        ← Pairwise alignments of mitogenomes vs reference
│   └── Genus_species_*_vs_ref.maf    ← LASTZ output: each mitogenome vs O. immaculatus ref
│
├── loci/                             ← All sequence data organized by locus and source
│   │
│   ├── ncbi_mito/                    ← Raw mito loci downloaded from NCBI GenBank
│   │   └── [locus]_ncbi.fa           ← 12S, 16S, COI, cytb, ND1-6, ATP6/8, COX2/3
│   │
│   ├── ncbi_nuclear/                 ← Raw nuclear loci downloaded from NCBI GenBank
│   │   └── [locus]_ncbi.fa           ← RAG1, MYH6, rhodopsin, ZIC1, etc.
│   │
│   ├── mito_extracted/               ← Mito loci extracted from YOUR assembled mitogenomes
│   │   └── [locus]_mito.fa           ← Same 15 loci, from mitogenome FASTAs
│   │
│   ├── ftol_mito/                    ← Mito loci sourced from Fish Tree of Life
│   │   └── 12s, 16s, coi, cytb, nd2, nd4 (subset available in FTOL)
│   │
│   ├── ftol_nuclear/                 ← Nuclear loci sourced from Fish Tree of Life
│   │   └── 22 loci (superset of ncbi_nuclear — some only in FTOL)
│   │
│   ├── merged/                       ← NCBI + extracted + FTOL combined per locus
│   │   ├── mito/[locus].fa           ← All sources merged, one file per mito locus
│   │   └── nuclear/[locus].fa        ← All sources merged, one file per nuclear locus
│   │
│   └── aligned/                      ← ★ Final aligned sequences ready for analysis
│       ├── mito/[locus]_aligned.fa   ← 15 mito loci aligned (MAFFT/MUSCLE)
│       └── nuclear/[locus]_aligned.fa ← 16 nuclear loci aligned
│
├── matrices/
│   └── morph_matrix                  ← Raw morphological character matrix
│
├── mitogenomes/                      ← Complete mitogenome FASTAs
│   ├── Genus_species_ACCESSION.fasta ← 13 mitogenomes used in analysis
│   ├── all_mitogenomes_aligned.fasta ← All mitogenomes aligned (full length)
│   ├── Ostracion_immaculatus_Ref_*   ← Reference mitogenome (used for LASTZ, extraction)
│   ├── additional_samples/           ← Extra mitogenomes not in main analysis
│   │   └── duplicates, old versions, extras
│   └── preliminary_files/            ← Early mitogenome alignment attempts
│       └── all_mitogenomes.fasta, .nex, total_evidence.nex
│
├── morpho/
│   └── Ostracioidea_K95_WT83.nex     ← ★ Morphological matrix used in total evidence analysis
│
├── r_cache/                          ← Cached R objects (speeds up re-running scripts)
│   ├── ftol_full.rds                 ← Cached full FTOL dataset
│   └── ncbi_results.rds              ← Cached NCBI query results
│
├── reference/                        ← Reference genome for alignment/annotation
│   ├── Ostracion_immaculatus_Ref_mitogenome_NC_009865.1.fa   ← Reference FASTA
│   └── Ostracion_immaculatus_Ref_mitogenome_NC_009865.1.gb   ← GenBank annotation file
│
├── results/                          ← ★ Figures and analysis scripts (presentable outputs)
│   ├── boxfish_data.qmd              ← R/Quarto script: data assembly & coverage analysis
│   ├── boxfishscript.qmd             ← R/Quarto script: main analysis & figures
│   ├── beast_tree_1.png              ← BEAST MCC tree figure
│   ├── iqtree_astral_tree_1.png      ← IQ-TREE / ASTRAL tree figure
│   └── data_coverage_by_source.png   ← Heatmap of data availability per species/locus
│
├── supermatrix/                      ← ★ Final concatenated matrices fed into BEAST
│   ├── supermatrix_26taxa.nex        ← Combined molecular + morpho matrix (26 taxa)
│   ├── molecular_26taxa.nex          ← Molecular partition only
│   ├── morpho_26taxa.nex             ← Morphological partition only
│   ├── supermatrix.nex               ← Earlier version
│   └── partitions.txt                ← Partition scheme: which bp = which locus/model
│
└── trees/                            ← ★ All final tree outputs
    ├── beast_MCC.tre                 ← (copy) Bayesian time-calibrated tree
    ├── boxfish_mito_iqtree.treefile  ← IQ-TREE mitogenome-only ML tree
    ├── astral_species_tree_final.tre ← ASTRAL coalescent species tree
    ├── boxfish_total.con.tre         ← Total evidence consensus
    └── morpho_topology.tre           ← Morphology-only topology
