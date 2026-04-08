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
