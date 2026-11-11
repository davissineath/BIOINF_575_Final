Function breakdown:
1. Load  
   read_chip(path) -> chip table, rsIDs only                                   A  
   read_pharmgkb(path) -> annotation table, columns renamed, every line kept   A    
   read_clinical_variants(path) -> rsID, gene, level, one drug per row         C  
   read_genotypes(path) -> {rsID: genotype}                                    A  
2. Analysis  
   split_column(table, column) -> one row per value of "warfarin, acenocoumarol" C  
   merge_chip_pharmgkb(chip, pgx) -> merged table (task 1)                     A  
   filter_associations(table, category="efficacy") -> tasks 2 and 3            B  
   coverage_per_gene(pgx_efficacy, chip) -> coverage table (task 3)            B  
   strongest_evidence(map_table, clinical) -> evidence table (task 4)          C  
   parse_alleles("AA + AG") -> ["AA", "AG"]                                    C  
   carrier_status("GA", ["AG"]) -> "yes"                                       C  
   person_report(evidence, genotypes) -> report table (task 5)                 B  
3. Plots  
   plot_coverage(coverage, path) -> the bar plot (task 3)                      B  
5. Main  
   main(data_dir, out_dir) -> runs everything, saves the files                 C  
