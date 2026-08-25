#!/bin/bash
set -e

# Combine in planta + in vitro metadata into one file first
python3 -c "
import pandas as pd
df1 = pd.read_csv('in_planta_master_metadata.csv')
df2 = pd.read_csv('in_vitro_master_metadata.csv')
pd.concat([df1, df2], ignore_index=True).to_csv('combined_master_metadata.csv', index=False)
"

python3 03a_assign_strain_identity.py \
  --metadata combined_master_metadata.csv \
  --seqtab-r1 \
    341F2_799F1_run2_invitro_R1.csv \
    341F2_799F1_run3_inplanta_R1.csv \
    341F2_799F2_run2_invitro_R1.csv \
    341F2_799F2_run3_inplanta_R1.csv \
    341F3_799F1_run2_invitro_R1.csv \
    341F3_799F1_run3_inplanta_R1.csv \
    341F3_799F2_run2_inplanta_R1.csv \
    341F6_799F2_run2_inplanta_R1.csv \
  --seqtab-r2 \
    341F2_799F1_run2_invitro_R2.csv \
    341F2_799F1_run3_inplanta_R2.csv \
    341F2_799F2_run2_invitro_R2.csv \
    341F2_799F2_run3_inplanta_R2.csv \
    341F3_799F1_run2_invitro_R2.csv \
    341F3_799F1_run3_inplanta_R2.csv \
    341F3_799F2_run2_inplanta_R2.csv \
    341F6_799F2_run2_inplanta_R2.csv \
  --ref-r1 monoculture_reference/monoculture_reference_R1.csv \
  --ref-r2 monoculture_reference/monoculture_reference_R2.csv \
  --output-dir strain_assignment
