# PyProteoMass_INTR-CADS-GP
PyProteoMass is an automated Python pipeline developed to bridge the gap between raw biological sequence data and genomic annotation files by parsing genome-scale proteomes, computing exact protein molecular weights, and mapping them to gene locus tags.

Key Architecture & Methodology: 
Streamed Parsing: Utilizes Biopython (Bio.SeqIO) generator streaming for memory-efficient sequence parsing and mass calculation.   
Dynamic Ingestion: Employs Pandas-driven case-insensitive header auto-detection and delimiter fallbacks to handle heterogeneous NCBI formats (.faa, .tsv, .csv, .ptt).   
Biochemical Accuracy: Applies condensation stoichiometry with water loss correction (subtracting $18.015\text{ Da}$ per peptide bond) and filters non-standard/ambiguous amino acid characters.   
Vectorized Performance: Uses $O(1)$ hash table lookups via Pandas vectorized mapping for rapid data integration.   

Performance & Validation:
Empirical Results: Successfully processed all 4,300+ protein-coding genes of Escherichia coli K-12 MG1655 in under 3 seconds with 100% locus tag matching accuracy.   Cross-Domain Scalability: Demonstrated native execution across both prokaryotic (E. coli) and eukaryotic (S. cerevisiae) datasets without requiring code modifications.   

Primary Output: 
The system generates standardized, two-column CSV datasets (locus_tag, Molecular_Weight_Da) that are fully optimized for direct ingestion into downstream proteomics workflows, gel electrophoresis validation, and metabolic modeling.
