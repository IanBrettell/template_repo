# Bash script to run on EBI cluster:

# conda activate snakemake
# snakemake \
#   --jobs 5000 \
#   --latency-wait 100 \
#   --cluster-config pilot_paper/snmk/config/cluster.json \
#   --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -n {cluster.n} -M {cluster.memory} -o {cluster.output} -e {cluster.error}' \
#   --keep-going \
#   --rerun-incomplete \
#   --use-conda \
#   -s pilot_paper/snmk/snakefile.py \
#   -p

# Import functions and packages

from os.path import join
import pandas as pd

# Load config file and provide content as config object

configfile: "pilot_paper/snmk/config/config.yaml"

# Load samples to process

SAMPLES = pd.read_csv(config["samples_file"], comment="#", skip_blank_lines=True, index_col=0)

# Globals

ASSAYS = ["open_field", "novel_object"]
QUADRANTS = ["q1", "q2", "q3", "q4"]

# Rules

rule all:
    input:
        expand("split/{assay}/{sample}_{quadrant}_{assay}.mp4",
            assay = ASSAYS,
            sample = SAMPLES.index,
            quadrant = QUADRANTS)

rule split:
    input:
        join(config["input_dir"], "{sample}.avi")
    params:
        samples_file = config["samples_file"]
    output:
        join(config["dlc_project_path"], "videos/{sample}_{quadrant}_{assay}.mp4")
    conda:
        config["python_env"]
    script:
        config["split_script"]

rule analyse_videos:
    input:
        join(config["dlc_project_path"], "videos/{sample}_{quadrant}_{assay}.mp4")
    params:
        dlc_project_path = config["dlc_project_path"]
    output:
        full = join(config["dlc_project_path"], "videos/{sample}_{quadrant}_{assay}") + config["dlc_scorer_name"] + "_full.pickle",
        meta = join(config["dlc_project_path"], "videos/{sample}_{quadrant}_{assay}") + config["dlc_scorer_name"] + "_meta.pickle"
    conda:
        config["dlc-cpu_env"]
    script:
        config["dlc_analyze_script"]
