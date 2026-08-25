#!/bin/bash
set -e

# 341F2_799F1_run2  (covers: in_vitro plate 1)
Rscript 02a_export_seqtab.R dada2_invitro_out/341F2_799F1_run2/seqtab_nochim_R1.rds 341F2_799F1_run2_invitro_R1.csv
Rscript 02a_export_seqtab.R dada2_invitro_out/341F2_799F1_run2/seqtab_nochim_R2.rds 341F2_799F1_run2_invitro_R2.csv

# 341F2_799F1_run3  (covers: in_planta plate 5)
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F2_799F1_run3/seqtab_nochim_R1.rds 341F2_799F1_run3_inplanta_R1.csv
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F2_799F1_run3/seqtab_nochim_R2.rds 341F2_799F1_run3_inplanta_R2.csv

# 341F2_799F2_run2  (covers: in_vitro plate 3)
Rscript 02a_export_seqtab.R dada2_invitro_out/341F2_799F2_run2/seqtab_nochim_R1.rds 341F2_799F2_run2_invitro_R1.csv
Rscript 02a_export_seqtab.R dada2_invitro_out/341F2_799F2_run2/seqtab_nochim_R2.rds 341F2_799F2_run2_invitro_R2.csv

# 341F2_799F2_run3  (covers: in_planta plate 4a, in_planta plate 4b)
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F2_799F2_run3/seqtab_nochim_R1.rds 341F2_799F2_run3_inplanta_R1.csv
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F2_799F2_run3/seqtab_nochim_R2.rds 341F2_799F2_run3_inplanta_R2.csv

# 341F3_799F1_run2  (covers: in_vitro plate 2)
Rscript 02a_export_seqtab.R dada2_invitro_out/341F3_799F1_run2/seqtab_nochim_R1.rds 341F3_799F1_run2_invitro_R1.csv
Rscript 02a_export_seqtab.R dada2_invitro_out/341F3_799F1_run2/seqtab_nochim_R2.rds 341F3_799F1_run2_invitro_R2.csv

# 341F3_799F1_run3  (covers: in_planta plate 3)
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F3_799F1_run3/seqtab_nochim_R1.rds 341F3_799F1_run3_inplanta_R1.csv
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F3_799F1_run3/seqtab_nochim_R2.rds 341F3_799F1_run3_inplanta_R2.csv

# 341F3_799F2_run2  (covers: in_planta plate 1)
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F3_799F2_run2/seqtab_nochim_R1.rds 341F3_799F2_run2_inplanta_R1.csv
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F3_799F2_run2/seqtab_nochim_R2.rds 341F3_799F2_run2_inplanta_R2.csv

# 341F6_799F2_run2  (covers: in_planta plate 2)
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F6_799F2_run2/seqtab_nochim_R1.rds 341F6_799F2_run2_inplanta_R1.csv
Rscript 02a_export_seqtab.R dada2_inplanta_out/341F6_799F2_run2/seqtab_nochim_R2.rds 341F6_799F2_run2_inplanta_R2.csv
