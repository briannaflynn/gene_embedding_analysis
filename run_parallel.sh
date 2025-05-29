#!/bin/bash

CUDA_VISIBLE_DEVICES=0 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_0_bin0.fasta" "/home/ubuntu/output/pair_0_bin0" &
CUDA_VISIBLE_DEVICES=1 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_10_bin2.fasta" "/home/ubuntu/output/pair_10_bin2" &
CUDA_VISIBLE_DEVICES=2 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_11_bin2.fasta" "/home/ubuntu/output/pair_11_bin2" &
CUDA_VISIBLE_DEVICES=3 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_12_bin2.fasta" "/home/ubuntu/output/pair_12_bin2" &
CUDA_VISIBLE_DEVICES=4 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_13_bin2.fasta" "/home/ubuntu/output/pair_13_bin2" &
CUDA_VISIBLE_DEVICES=5 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_14_bin2.fasta" "/home/ubuntu/output/pair_14_bin2" &
CUDA_VISIBLE_DEVICES=6 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_15_bin3.fasta" "/home/ubuntu/output/pair_15_bin3" &
CUDA_VISIBLE_DEVICES=7 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_16_bin3.fasta" "/home/ubuntu/output/pair_16_bin3" &
CUDA_VISIBLE_DEVICES=0 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_17_bin3.fasta" "/home/ubuntu/output/pair_17_bin3" &
CUDA_VISIBLE_DEVICES=1 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_18_bin3.fasta" "/home/ubuntu/output/pair_18_bin3" &
CUDA_VISIBLE_DEVICES=2 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_19_bin3.fasta" "/home/ubuntu/output/pair_19_bin3" &
CUDA_VISIBLE_DEVICES=3 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_1_bin0.fasta" "/home/ubuntu/output/pair_1_bin0" &
CUDA_VISIBLE_DEVICES=4 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_20_bin4.fasta" "/home/ubuntu/output/pair_20_bin4" &
CUDA_VISIBLE_DEVICES=5 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_21_bin4.fasta" "/home/ubuntu/output/pair_21_bin4" &
CUDA_VISIBLE_DEVICES=6 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_22_bin4.fasta" "/home/ubuntu/output/pair_22_bin4" &
CUDA_VISIBLE_DEVICES=7 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_23_bin4.fasta" "/home/ubuntu/output/pair_23_bin4" &
CUDA_VISIBLE_DEVICES=0 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_24_bin4.fasta" "/home/ubuntu/output/pair_24_bin4" &
CUDA_VISIBLE_DEVICES=1 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_2_bin0.fasta" "/home/ubuntu/output/pair_2_bin0" &
CUDA_VISIBLE_DEVICES=2 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_3_bin0.fasta" "/home/ubuntu/output/pair_3_bin0" &
CUDA_VISIBLE_DEVICES=3 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_4_bin0.fasta" "/home/ubuntu/output/pair_4_bin0" &
CUDA_VISIBLE_DEVICES=4 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_5_bin1.fasta" "/home/ubuntu/output/pair_5_bin1" &
CUDA_VISIBLE_DEVICES=5 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_6_bin1.fasta" "/home/ubuntu/output/pair_6_bin1" &
CUDA_VISIBLE_DEVICES=6 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_7_bin1.fasta" "/home/ubuntu/output/pair_7_bin1" &
CUDA_VISIBLE_DEVICES=7 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_8_bin1.fasta" "/home/ubuntu/output/pair_8_bin1" &
CUDA_VISIBLE_DEVICES=0 colabfold_batch --model-type alphafold2_multimer_v3 "alphafold2multimer_inputs/pair_9_bin1.fasta" "/home/ubuntu/output/pair_9_bin1" &
wait
