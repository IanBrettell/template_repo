#!/usr/bin/env python3

# Import libaries

import os
os.environ["DLClight"]="True" # suppress GUI

import deeplabcut

# Setup project variables

config_path = os.path.join(snakemake.params.dlc_project_path, 'config.yaml')
# dlc_project_path = 'dlc/pp_proj-ian_brettell-2020-11-17'
# config_path = os.path.join(dlc_project_path, 'config.yaml')

target_video = snakemake.input[0]
# target_video = 'dlc/pp_proj-ian_brettell-2020-11-17/videos/20190611_1410_icab_kaga_R_q4_open_field.mp4'

# Analyse video

deeplabcut.analyze_videos(config_path, target_video, videotype='mp4', save_as_csv=True)
