library(ape)
library(phangorn)
library(tidyverse)
library(ggtree)

# ══════════════════════════════════════════════════════════════════════
# 0. MASTER LOOKUP TABLES (defined once)
# ══════════════════════════════════════════════════════════════════════

taxonomy <- data.frame(
  species_full = c(
    "Anoplocapros amygdaloides", "Anoplocapros lenticularis", "Anoplocapros inermis",
    "Strophiurichthys robustus", "Aracana aurita", "Aracana ornata",
    "Capropygia unistriata", "Kentrocapros aculeatus", "Polyplacapros tyleri",
    "Acanthostracion guineensis", "Acanthostracion notacanthus",
    "Acanthostracion polygonius", "Acanthostracion quadricornis",
    "Lactophrys bicaudalis", "Lactophrys trigonus", "Lactophrys triqueter",
    "Lactoria cornuta", "Lactoria diaphana", "Lactoria fornasini", "Lactoria paschae",
    "Ostracion cubicus", "Ostracion cyanurus", "Ostracion immaculatus",
    "Ostracion meleagris", "Ostracion nasus", "Rhynchostracion rhinorhynchus",
    "Ostracion solorensis", "Ostracion whitleyi",
    "Tetrosomus concatenatus", "Tetrosomus gibbosus", "Ostracion reipublicae",
    "Triodon macropterus"
  ),
  morpho_code = c(
    NA, "Anoplocapros", NA,
    "Strophiurichthys", "Aracana", NA,
    "Capropygia", "Kentrocapros", "Polyplacapros",
    "acant_guin", "acant_nota", "acant_poly", "acant_quad",
    "lactp_bica", "lactp_trig", "lactp_triq",
    "lacto_corn", "lacto_diap", "lacto_forn", NA,
    "ostr_cubic", "ostr_cyanu", NA, "ostr_melea",
    "ostr_nasus", "Ostr_rhino", NA, "ostr_whitl",
    "tetr_conca", "tetr_gibbo", "tetr_repub",
    "Triodon"
  ),
  mitogenome_accession = c(
    NA, "NC_011319.1", NA,
    NA, NA, NA, NA, "NC_009864.1", NA, NA, NA,
    "PV662072.1", "PV662073.1",
    "OP057014.1", "OP035078.2", NA,
    "NC_057981.1", "NC_011330.1", NA, NA,
    "NC_072284.1", NA, "NC_009865.1", "OP035246.1",
    NA, "NC_030005.1", NA, NA, NA, NA, NA,
    "NC_009859.1"
  ),
  ftol_name = c(
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    "Acanthostracion_polygonius", "Acanthostracion_quadricornis",
    "Lactophrys_bicaudalis", "Lactophrys_trigonus", "Rhinesomus_triqueter",
    "Lactoria_cornuta", "Lactoria_diaphana", "Tetrosomus_fornasini", NA,
    "Ostracion_cubicus", NA, "Ostracion_immaculatus", "Ostracion_meleagris",
    "Rhynchostracion_nasus", "Ostracion_rhinorhynchos",
    "Ostracion_solorensis", "Ostracion_whitleyi",
    "Tetrosomus_concatenatus", "Tetrosomus_gibbosus", "Tetrosomus_reipublicae",
    "Triodon_macropterus"
  ),
  family = c(rep("Aracanidae", 9), rep("Ostraciidae", 22), "Triodontidae"),
  stringsAsFactors = FALSE
) |>
  mutate(species_underscore = gsub(" ", "_", species_full))

# Convenience named vectors for fast tip relabelling
acc_to_species  <- setNames(taxonomy$species_full,  taxonomy$mitogenome_accession)
code_to_species <- setNames(taxonomy$species_full,  taxonomy$morpho_code)
acc_to_code     <- setNames(taxonomy$morpho_code,   taxonomy$mitogenome_accession)


# ══════════════════════════════════════════════════════════════════════
# 1. READ DATA (each file read once)
# ══════════════════════════════════════════════════════════════════════

morpho_raw <- read.nexus.data("morpho/Ostracioidea_K95_WT83.nex")
morpho_mat <- do.call(rbind, morpho_raw)

dna_data <- read.dna(
  "muscle_results/muscle-I20260206-203642-0302-71998231-p1m.sqr.fa",
  format = "fasta"
)

mol_aln <- read.dna(
  "boxfish-seqs/mitogenomes/all_mitogenomes_aligned.fasta",
  format = "fasta"
)
rownames(mol_aln) <- gsub(" .*", "", rownames(mol_aln))   # strip to accession

# Pre-computed: Ostraciidae nuclear alignment already has Triodon appended
# (saves re-reading the full ~24k-col FToL alignment just to extract one row)
combined_aln <- read.dna(
  "~/Downloads/Boxfish/ftol/Ostraciidae_with_triodon.phylip",
  format = "sequential"
)

# ostrax is still needed for coverage bookkeeping
ostrax <- combined_aln[rownames(combined_aln) != "Triodon_macropterus", ]


# ══════════════════════════════════════════════════════════════════════
# 2. MERGE DNA + LOOKUP
# ══════════════════════════════════════════════════════════════════════

dna_named <- as.data.frame(as.character(dna_data)) |>
  rownames_to_column("ID") |>
  left_join(taxonomy[, c("mitogenome_accession", "morpho_code")],
            by = c("ID" = "mitogenome_accession")) |>
  filter(!is.na(morpho_code)) |>
  column_to_rownames("morpho_code") |>
  select(-ID)

unmatched <- rownames(as.character(dna_data))[!rownames(as.character(dna_data)) %in% taxonomy$mitogenome_accession]
if (length(unmatched)) message("Unmatched IDs: ", paste(unmatched, collapse = ", "))

shared_taxa    <- intersect(rownames(dna_named), rownames(morpho_mat))
message(length(shared_taxa), " shared taxa found")
combined_mat   <- cbind(morpho_mat[shared_taxa, ], dna_named[shared_taxa, ])


# ══════════════════════════════════════════════════════════════════════
# 3. SPLIT ALIGNMENT INTO LOCI
# ══════════════════════════════════════════════════════════════════════

aln <- read.dna("Ostraciidae.phylip", format = "sequential")

locus_ranges <- list(
  # Mitochondrial
  s12s      = 1:979,
  s16s      = 980:1756,
  coi       = 2292:2973,
  cytb      = 2974:4114,
  nd2       = 8822:9868,
  nd4       = 9869:11401,
  # Nuclear
  s4c4      = 1757:2291,
  enc1      = 4115:4954,
  ficd      = 4955:5680,
  glyt      = 5681:6613,
  hoxc6a    = 6614:7243,
  kiaa1239  = 7244:7987,
  myh6      = 7988:8821,
  panx2     = 11402:12141,
  plagl2    = 12142:12951,
  ptr       = 12952:13662,
  rag1      = 13663:15116,
  rag2      = 15117:16340,
  rhodopsin = 16341:17264,
  ripk4     = 17265:17909,
  sh3px3    = 17910:18632,
  sidkey    = 18633:19727,
  sreb2     = 19728:20714,
  svep1     = 20715:21521,
  tbr1      = 21522:22349,
  vcpip     = 22350:23114,
  zic1      = 23115:24143
)

mito_loci    <- c("s12s","s16s","coi","cytb","nd2","nd4")
nuclear_loci <- setdiff(names(locus_ranges), mito_loci)

dir.create("loci/mito",    recursive = TRUE, showWarnings = FALSE)
dir.create("loci/nuclear", recursive = TRUE, showWarnings = FALSE)

for (nm in names(locus_ranges)) {
  subdir <- if (nm %in% mito_loci) "loci/mito" else "loci/nuclear"
  fname  <- file.path(subdir, paste0(nm, ".fa"))
  write.dna(aln[, locus_ranges[[nm]]], fname, format = "fasta")
}

message("Done! ", ncol(aln), " total columns split into ", length(locus_ranges),
        " loci across loci/mito/ and loci/nuclear/")


# ══════════════════════════════════════════════════════════════════════
# 4. MITOGENOME TREE — read, root, relabel, plot
# ══════════════════════════════════════════════════════════════════════

relabel_tips <- function(tree, lookup_vec) {
  tree$tip.label <- lookup_vec[tree$tip.label]
  tree
}

mito_tree <- read.tree("~/Downloads/Boxfish/boxfish_mito_iqtree.treefile") |>
  relabel_tips(acc_to_species) |>
  (\(t) root(t, grep("Triodon", t$tip.label, value = TRUE), resolve.root = TRUE))()

par(mar = c(1, 1, 2, 1))
plot(mito_tree,
     main = "Boxfish Mitogenome Phylogeny (GTR+G, 1000 bootstrap)",
     cex = 0.8, label.offset = 0.005)
nodelabels(mito_tree$node.label,
           cex = 0.6, frame = "none", adj = c(1.2, -0.3), col = "blue")
add.scale.bar()


# ══════════════════════════════════════════════════════════════════════
# 5. MORPHOLOGICAL PARSIMONY TREE
# ══════════════════════════════════════════════════════════════════════

morpho_phydat <- phyDat(morpho_mat, type = "USER", levels = c("0","1","2","3"))

set.seed(42)
morpho_mp <- pratchet(morpho_phydat, trace = 0, minit = 100, k = 10) |>
  (\(t) root(t, "Aracana", resolve.root = TRUE))()


# ══════════════════════════════════════════════════════════════════════
# 6. PRUNE & COMPARE TREES (shared taxa)
# ══════════════════════════════════════════════════════════════════════

# Map mito accessions → morpho codes for the mito tree
mito_tree_coded <- read.tree("~/Downloads/Boxfish/boxfish_mito_iqtree.treefile")
mito_tree_coded$tip.label <- acc_to_code[mito_tree_coded$tip.label]

shared <- intersect(mito_tree_coded$tip.label, names(morpho_raw))
cat("Shared taxa:", length(shared), "\n")

mito_pruned  <- keep.tip(mito_tree_coded, shared)
morpho_pruned <- keep.tip(morpho_mp, shared)

par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot(morpho_pruned, main = "Morphological Parsimony\n(shared taxa)", cex = 0.8)
plot(mito_pruned,   main = "Mitogenome ML\n(shared taxa)",            cex = 0.8)

# Save morpho tree with accession IDs for MrBayes
morpho_renamed <- morpho_pruned
morpho_renamed$tip.label <- taxonomy$mitogenome_accession[
  match(morpho_renamed$tip.label, taxonomy$morpho_code)
]
write.tree(morpho_renamed, "morpho_topology.tre")


# ══════════════════════════════════════════════════════════════════════
# 7. ADD TRIODON TO NUCLEAR ALIGNMENT
# ══════════════════════════════════════════════════════════════════════

triodon_idx  <- grep("Triodon_macropterus", rownames(ftol_full))
combined_aln <- rbind(ostrax, ftol_full[triodon_idx, ])
rownames(combined_aln)[nrow(combined_aln)] <- "Triodon_macropterus"

write.dna(combined_aln,
          file   = "~/Downloads/Boxfish/ftol/Ostraciidae_with_triodon.phylip",
          format = "sequential",
          colsep = "", nbcol = -1)

cat("Taxa:", nrow(combined_aln), " | Characters:", ncol(combined_aln), "\n")


# ══════════════════════════════════════════════════════════════════════
# 8. BUILD TOTAL-EVIDENCE NEXUS FOR MRBAYES
# ══════════════════════════════════════════════════════════════════════

all_taxa     <- taxonomy$morpho_code
n_mol_chars  <- ncol(mol_aln)
n_morph_chars <- ncol(morpho_mat)

# Molecular matrix — fill missing with "?"
mol_combined <- matrix("?", nrow = nrow(taxonomy), ncol = n_mol_chars,
                       dimnames = list(all_taxa, NULL))
for (i in seq_len(nrow(taxonomy))) {
  acc <- taxonomy$mitogenome_accession[i]
  if (!is.na(acc) && acc %in% rownames(mol_aln))
    mol_combined[i, ] <- as.character(mol_aln[acc, ])
}

# Morpho matrix — fill missing with "?"
morpho_combined <- matrix("?", nrow = nrow(taxonomy), ncol = n_morph_chars,
                           dimnames = list(all_taxa, NULL))
for (i in seq_len(nrow(taxonomy))) {
  code <- taxonomy$morpho_code[i]
  if (!is.na(code) && code %in% names(morpho_raw))
    morpho_combined[i, ] <- morpho_raw[[code]]
}

n_total     <- n_mol_chars + n_morph_chars
outgrp_name <- all_taxa[grep("riodon", all_taxa, ignore.case = TRUE)]

con <- file("boxfish-seqs/mitogenomes/total_evidence.nex", "w")
writeLines(c(
  "#NEXUS", "",
  "begin data;",
  paste0("  dimensions ntax=", nrow(taxonomy), " nchar=", n_total, ";"),
  paste0("  format datatype=mixed(dna:1-", n_mol_chars,
         ", standard:", n_mol_chars + 1, "-", n_total,
         ") interleave=no missing=? gap=-;"),
  "  matrix"
), con)

for (i in seq_len(nrow(taxonomy))) {
  writeLines(paste0("  ", all_taxa[i], "  ",
                    paste(mol_combined[i,],    collapse = ""),
                    paste(morpho_combined[i,], collapse = "")), con)
}

writeLines(c(
  "  ;", "end;", "",
  "begin mrbayes;",
  "  set autoclose=yes nowarn=yes;", "",
  paste0("  charset molecular = 1-", n_mol_chars, ";"),
  paste0("  charset morpho = ", n_mol_chars + 1, "-", n_total, ";"),
  "  partition totalevidence = 2: molecular, morpho;",
  "  set partition = totalevidence;", "",
  "  lset applyto=(1) nst=6 rates=gamma ngammacat=4;",
  "  lset applyto=(2) nst=1 rates=gamma ngammacat=4 coding=variable;", "",
  "  prset applyto=(1) statefreqpr=dirichlet(1,1,1,1);",
  "  prset applyto=(2) symdirihyperpr=fixed(infinity);",
  "  prset applyto=(all) ratepr=variable;", "",
  "  unlink statefreq=(all) revmat=(all) shape=(all);", "",
  paste0("  outgroup ", outgrp_name, ";"),
  "  mcmc ngen=2000000 nchains=4 nruns=2",
  "       samplefreq=500 printfreq=500",
  "       savebrlens=yes filename=boxfish_total;",
  "  sump burnin=500;",
  "  sumt burnin=500;",
  "end;"
), con)
close(con)

cat("Total-evidence nexus written.\n")


# ══════════════════════════════════════════════════════════════════════
# 8. THREE-TREE COMPARISON PLOT
# ══════════════════════════════════════════════════════════════════════

nuclear_tree <- read.tree("~/Downloads/Boxfish/ftol/nuclear_iqtree.treefile") |>
  (\(t) root(t, "Triodon_macropterus", resolve.root = TRUE))()
nuclear_tree$tip.label <- gsub("_", " ", nuclear_tree$tip.label)

morpho_tree <- read.tree("~/Downloads/Boxfish/morpho_topology.tre") |>
  relabel_tips(acc_to_species) |>
  (\(t) root(t, "Kentrocapros aculeatus", resolve.root = TRUE))()

par(mfrow = c(1, 3))
plot(mito_tree,    main = "Mitogenome (IQ-TREE)",        cex = 0.7)
nodelabels(mito_tree$node.label,    cex = 0.5, frame = "none", adj = c(1.2, -0.3))
plot(nuclear_tree, main = "Nuclear loci (IQ-TREE)",      cex = 0.7)
nodelabels(nuclear_tree$node.label, cex = 0.5, frame = "none", adj = c(1.2, -0.3))
plot(morpho_tree,  main = "Morphology (Parsimony)",      cex = 0.7)


# ══════════════════════════════════════════════════════════════════════
# 10. DATA COVERAGE HEATMAP
# ══════════════════════════════════════════════════════════════════════

mito_species    <- na.omit(taxonomy$species_full[!is.na(taxonomy$mitogenome_accession)])
ftol_species    <- gsub("_", " ", na.omit(taxonomy$ftol_name))
morpho_species  <- na.omit(taxonomy$species_full[!is.na(taxonomy$morpho_code)])
target_species  <- taxonomy$species_full

coverage <- taxonomy |>
  transmute(
    Species      = species_full,
    Mitogenome   = species_full %in% mito_species,
    FToL_nuclear = species_full %in% ftol_species,
    Morphology   = species_full %in% morpho_species
  )

coverage |>
  pivot_longer(-Species, names_to = "Data_type", values_to = "Present") |>
  mutate(
    Data_type = factor(Data_type, levels = c("Mitogenome","FToL_nuclear","Morphology")),
    Species   = factor(Species, levels = rev(target_species))
  ) |>
  ggplot(aes(x = Data_type, y = Species, fill = Present)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#2d6a4f", "FALSE" = "#f0f0f0"),
                    labels  = c("TRUE" = "Present", "FALSE" = "Absent")) +
  labs(title = "Data coverage across taxa", x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic", size = 8))


# ══════════════════════════════════════════════════════════════════════
# 9. NCBI COVERAGE HEATMAP (load pre-computed results)
# ══════════════════════════════════════════════════════════════════════

# Previously saved by the rentrez search — no need to re-run
coverage_ncbi <- read.csv("~/Downloads/Boxfish/ncbi_coverage.csv",
                           stringsAsFactors = FALSE)

locus_order <- c(
  "COI","cytb","16S","12S","ND2","ND4","ND1","ND3","ND4L","ND5","ND6",
  "COX2","COX3","ATP6","ATP8",
  "4c4","enc1","ficd","glyt","hoxc6a","kiaa1239","myh6","panx2","plagl2",
  "ptr","rag1","rag2","rhodopsin","ripk4","sh3px3","sidkey","sreb2",
  "svep1","tbr1","vcpip","zic1"
)

coverage_ncbi |>
  mutate(
    locus   = factor(locus,   levels = locus_order),
    species = factor(species, levels = rev(target_species))
  ) |>
  ggplot(aes(x = locus, y = species, fill = found)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#2d6a4f", "FALSE" = "#f0f0f0"),
                    labels  = c("TRUE" = "Present", "FALSE" = "Absent")) +
  labs(title = "NCBI data availability across taxa and loci",
       x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(face = "italic", size = 8))
