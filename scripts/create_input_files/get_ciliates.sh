# get to scripts
cd /storage/home/hcoda1/5/cbeuscher3/r-alind6-0/orthologous_euks/SM/scripts/create_input_files
INPUT_DIR=/storage/home/hcoda1/5/cbeuscher3/r-alind6-0/orthologous_euks/SM/input_files


wget -O RCG_450.gtf.zip "https://ndownloader.figshare.com/files/60734557"
wget -O RCG_450.pep.zip "https://ndownloader.figshare.com/files/60734557"


#insert code to download sample key


    python3 filter_genomes.py \\
        "$INPUT_DIR/RCG_450_summary.xlsx" \\
        "$INPUT_DIR/rep_RCG_450_summary.xlsx"

#generate sample ids from sample key
# ran sucessfully on 2024-06-05
    python3 generate_sample_ids.py \\
        "$INPUT_DIR/RCG_450_summary.xlsx" \\
        "$INPUT_DIR/sample_key.tsv" \\
        "$INPUT_DIR/sample_species_key.tsv" \\ #sample_id /t species
        "$INPUT_DIR/genomename_sampleid_key.tsv" # og_name /t sample_id

# append to exisiting sample key
cat "$INPUT_DIR/sample_species_key.tsv" >> sample_key.tsv

#rename genomes in sample key as new sample ids
    python3 rename_xlsx_genomes.py \\
        "$INPUT_DIR/RCG_450_summary.xlsx" \\
        "$INPUT_DIR/genomename_sampleid_key.tsv" \\
        "$INPUT_DIR/RCG_450_summary_renamed.xlsx"

# Parse taxonomy from sample key column in same format as grouping key
    python3 parse_genome_xlsx.py \\
        "$INPUT_DIR/RCG_450_summary.xlsx" \\
        "$INPUT_DIR/new_groupings.tsv"

# preserve original copy of grouping key

cp "$INPUT_DIR/grouping_key.tsv" "$INPUT_DIR/og_grouping_key.tsv"
# add output to grouping key
    python3 combine_taxonomy.py \\
        "$INPUT_DIR/og_grouping_key.tsv" \\
        "$INPUT_DIR/new_groupings.tsv" \\
        "$INPUT_DIR/grouping_key.tsv"


    python3 combine_taxonomy.py \\
        "$INPUT_DIR/og_grouping_key.tsv" \\
        "$INPUT_DIR/new_groupings.tsv" \\
        "$INPUT_DIR/input_files/HA.txt" \\
        "$INPUT_DIR//aerobic.txt" \\
        "$INPUT_DIR/grouping_key.tsv" 


python combine_taxonomy.py ../../../input_files/before_rumen/og_grouping_key.tsv ../../../input_files/new_groupings.tsv ../../../input_files/HA.txt ../../../input_files/aerobic.txt ../../../input_files/grouping_key.tsv




#filter for representative genomes
    python3 filter_genomes.py \\
        "$INPUT_DIR/RCG_450_summary.xlsx" \\
        "$INPUT_DIR/rep_RCG_450_summary.xlsx"

#generate sample ids from sample key
    python3 generate_sample_ids.py \\
        "$INPUT_DIR/rep_RCG_450_summary.xlsx" \\
        "$INPUT_DIR/sample_key.tsv" \\
        "$INPUT_DIR/sample_species_key.tsv" \\ #sample_id /t species
        "$INPUT_DIR/genomename_sampleid_key.tsv" # og_name /t sample_id

# append to exisiting sample key
cat "$INPUT_DIR/sample_species_key.tsv" >> "$INPUT_DIR/sample_key.tsv"

#rename genomes in sample key as new sample ids
    python3 rename_xlsx_genomes.py \\
        "$INPUT_DIR/rep_RCG_450_summary.xlsx" \\
        "$INPUT_DIR/genome_sampleid_key.tsv.tsv" \\
        "$INPUT_DIR/new_RCG_450_summary.xlsx"

# Parse taxonomy from sample key column in same format as grouping key
    python3 parse_genome_xlsx.py \\
        "$INPUT_DIR/new_RCG_450_summary.xlsx" \\
        "$INPUT_DIR/new_groupings.tsv"

# preserve original copy of grouping key

cp "$INPUT_DIR/grouping_key.tsv" "$INPUT_DIR/og_grouping_key.tsv"
# add output to grouping key
    python3 combine_taxonomy.py \\
        "$INPUT_DIR/og_grouping_key.tsv" \\
        "$INPUT_DIR/new_groupings.tsv" \\
        "$INPUT_DIR/grouping_key.tsv"

    python3 combine_taxonomy.py \\
        "$INPUT_DIR/og_grouping_key.tsv" \\
        "$INPUT_DIR/new_groupings.tsv" \\
        "$INPUT_DIR/input_files/HA.txt" \\
        "$INPUT_DIR//aerobic.txt" \\
        "$INPUT_DIR/grouping_key.tsv"


            python3 rename_and_filter_faa.py \\
        "$INPUT_DIR/genome_sampleid_key.tsv" \\
        "$SOURCE_DIR" \\
        "$PRO_DIR"  



        snakemake -s Snakefile.smk --unlock \
    --config genome_dir=/storage/home/hcoda1/5/cbeuscher3/r-alind6-0/orthologous_euks/SM/fastas/fna protein_dir=/storage/home/hcoda1/5/cbeuscher3/r-alind6-0/orthologous_euks/SM/fastas/faa

