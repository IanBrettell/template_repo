# import functions and packages
from os.path import join
import pandas as pd

# load config file and provide content as config object
configfile: "config/config.yaml"

# load all samples we want to process
samples = pd.read_csv (config["sample_file"], comment="#", skip_blank_lines=True, index_col=0)

# Globals
ASSAYS = ["open_field", "novel_object"]
QUADRANTS = ["q1", "q2", "q3", "q4"]

# Rules
rule all:
    input:
        expand("split/{sample}_{quadrant}.avi",
            sample = samples.index,
            quadrant = QUADRANTS)

rule split:
    input:
        join(config["input_dir"], "{sample}.avi")
    output:
        "split/{sample}_{quadrant}.avi"
    shell:
        pilot_paper/code/scripts/20201111_split_videos.py \
            --in_file {input} \
  --start 2889\
  --end 22790\
  --quadrant q1 \
  --out_dir tmp
