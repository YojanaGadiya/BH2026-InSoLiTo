# Merged software pairs

- pipelines merged: 158
- distinct directed valid pairs: 14002 (12277 without utilities/MultiQC)
- valid pairs found in 2+ pipelines: 1324
- valid pairs by direction_type: {'one_direction': 6655, 'both_directions': 1124, 'parallel': 6223}
- pairs by best pair_type: {'data_flow': 7520, 'parallel_co_use': 6223, 'gate_metadata': 259, 'rejected': 6}

## pair_type

- `data_flow` (rank 1)
- `gate_metadata` (rank 2)
- `parallel_co_use` (rank 3)
- `rejected` (rank 4)

## Most frequent valid pairs (analysis tools only)

| source | target | pipelines | process edges | direction | best pair_type |
|---|---|---|---|---|---|
| bowtie2 | samtools | 25 | 116 | both_directions | data_flow |
| bwa | samtools | 22 | 77 | both_directions | data_flow |
| fastp | fastqc | 21 | 36 | one_direction | data_flow |
| samtools | picard | 20 | 265 | both_directions | data_flow |
| samtools | bedtools | 18 | 116 | both_directions | data_flow |
| samtools | fastqc | 18 | 9 | one_direction | data_flow |
| bowtie2 | fastqc | 16 | 5 | one_direction | data_flow |
| picard | samtools | 16 | 154 | both_directions | data_flow |
| bwa | fastqc | 15 | 0 | parallel | parallel_co_use |
| fastqc | bwa | 15 | 0 | parallel | parallel_co_use |
| minimap2 | samtools | 14 | 136 | both_directions | data_flow |
| trim-galore | fastqc | 14 | 7 | one_direction | data_flow |
| umi_tools | samtools | 14 | 244 | both_directions | data_flow |
| samtools | bcftools | 13 | 231 | both_directions | data_flow |
| fastqc | bowtie2 | 13 | 0 | parallel | parallel_co_use |
| fastp | bowtie2 | 12 | 124 | one_direction | data_flow |
| fastp | samtools | 12 | 151 | both_directions | data_flow |
| samtools | gatk4 | 12 | 342 | both_directions | data_flow |
| samtools | umi_tools | 12 | 118 | both_directions | data_flow |
| fastqc | samtools | 12 | 0 | parallel | parallel_co_use |
| fastp | bwa-mem2 | 11 | 84 | one_direction | data_flow |
| gffread | samtools | 11 | 31 | one_direction | data_flow |
| samtools | bwa | 11 | 55 | both_directions | data_flow |
| samtools | mosdepth | 11 | 92 | one_direction | data_flow |
| fastqc | trim-galore | 11 | 0 | parallel | parallel_co_use |
| bbmap | fastqc | 10 | 4 | one_direction | data_flow |
| bwa-mem2 | samtools | 10 | 79 | both_directions | data_flow |
| cutadapt | fastqc | 9 | 9 | one_direction | data_flow |
| fastp | bbmap | 9 | 22 | both_directions | data_flow |
| fastp | bwa | 9 | 70 | one_direction | data_flow |
| gatk4 | bcftools | 9 | 52 | both_directions | data_flow |
| samtools | subread | 9 | 46 | both_directions | data_flow |
| samplesheet_check | fastqc | 9 | 10 | one_direction | data_flow |
| bbmap | samtools | 8 | 65 | both_directions | data_flow |
| fastp | kraken2 | 8 | 20 | one_direction | data_flow |
| star | samtools | 8 | 12 | both_directions | data_flow |
| trim-galore | bowtie2 | 8 | 26 | both_directions | data_flow |
| trim-galore | bwa | 8 | 13 | one_direction | data_flow |
| fastqc | fastp | 8 | 0 | parallel | parallel_co_use |
| fastqc | minimap2 | 8 | 0 | parallel | parallel_co_use |
