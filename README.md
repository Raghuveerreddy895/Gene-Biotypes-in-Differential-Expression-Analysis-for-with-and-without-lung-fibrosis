# Gene Biotypes in Differential Expression Analysis for with and without lung fibrosis
n this project, we analyze genes based on both their expression behavior (Upregulated, Downregulated, Unexpressed, or Not Significant) and their functional biotype as annotated in Ensembl.  This helps us not only identify which genes are changing between conditions, but also understand what kinds of genes they are — whether they encode proteins, regulate other genes, or belong to small RNA classes with specialized roles.
# Gene Biotypes in DEG Analysis

Data sets accesion ID GSE213001

This project classifies genes based on both **differential expression status** 
(Upregulated, Downregulated, Unexpressed, Not Significant) 
and **Ensembl biotype**.

## Biotype Definitions

- **protein_coding** → Encodes proteins; typically main drivers of cellular function.  
- **lncRNA (long non-coding RNA)** → Non-protein coding transcripts (>200nt); 
  regulate transcription, splicing, chromatin remodeling, and can act as oncogenes/tumor suppressors.  
- **antisense** → lncRNA transcribed opposite to a coding gene; regulates stability and translation.  
- **processed_transcript** → Non-coding RNA without strong coding potential.  
- **sense_intronic** → lncRNA located within introns of coding genes.  
- **sense_overlapping** → lncRNA overlapping a coding gene in the same strand.  
- **3prime_overlapping_ncRNA** → Rare non-coding transcripts overlapping 3′ ends of protein-coding genes.  
- **macro_lncRNA** → Very long lncRNAs with regulatory functions.  
- **bidirectional_promoter_lncRNA** → lncRNAs transcribed in the opposite direction of a protein-coding gene promoter.  
- **pseudogene** → Resemble protein-coding genes but are non-functional due to mutations.  
- **snoRNA (small nucleolar RNA)** → Guide chemical modifications of rRNA.  
- **miRNA (microRNA)** → Short RNAs regulating mRNA degradation/translation.  
- **snRNA (small nuclear RNA)** → Components of the spliceosome for RNA splicing.  
- **rRNA (ribosomal RNA)** → Structural and functional components of ribosomes.  

## Key Files

- `DEGs_up_down_unexpressed_notsig.csv` → All genes with status & biotype.  
- `lncRNA_up_down_unexpressed_notsig.csv` → Long non-coding RNAs (lncRNAs) only.  
- `lncRNA_upregulated.csv` → Upregulated lncRNAs.  
- `lncRNA_downregulated.csv` → Downregulated lncRNAs.  
- `lncRNA_unexpressed.csv` → lncRNAs with zero counts.  
- `lncRNA_not_significant.csv` → lncRNAs with no significant differential expression.  
