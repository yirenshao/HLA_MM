# HLA_MM

## Description: 
The pipeline takes as input Donor and Host HLA typing (4 digits), HLA proteome - looks for peptides derived from HLA major mismatch and checks if peptides are present in remaining HLA Proteomes, if not present then predicts binding stability using HLAthena1.0 with the Donor HLAs (Peptide-HLA binding prediction), HLAthena1.0 post processing: picks out peptides that have a binding stability of <0.5. 
