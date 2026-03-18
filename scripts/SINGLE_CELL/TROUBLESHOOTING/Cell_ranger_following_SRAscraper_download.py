#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:38:36 2025

@author: jlehle
"""
#%% Read in the SRAscraper config.yaml to get the working dir
# This section is the only thing that need to be automated to bring it into the snakmake pipeline
import os
import yaml
#Change these two lines when you make up the CLI
os.chdir('/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer')
with open('config_NMF_v0.1.1.yaml', 'r') as config:
    config_yaml = yaml.safe_load(config)

#%% Go load the meatdata of all the files downloaded
working_dir = config_yaml['output_dir']
import pickle

with open(os.path.join(working_dir, 'metadata/dictionary_file.pkl'), 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

#%% Yehaw, cellranger time

import subprocess
from os import cpu_count
import pandas as pd
import traceback

NCPUS = cpu_count()
sample=("_S1_L001_")


# List of all possible chemistry options to try (in order of likelihood)
chemistry_options = [
    "auto",  # Try auto first
    "SC3Pv4",
    "SC3Pv3",
    "SC3Pv2",
    "SC3Pv3HT",
    "SC3Pv3LT",
    "threeprime",
    "fiveprime",
    "SC5P-PE-v3",
    "SC5P-PE",
    "SC5P-R2-v3",
    "SC5P-R2",
    "SC3Pv1",  # These last two can't be auto-detected
    "ARC-v1"
]

successful_runs = {key: pd.DataFrame(columns=gse_dict[key].columns) for key in gse_dict.keys()}

def format_error_output(result):
    """Format detailed error output from subprocess result"""
    error_lines = []
    if result.stdout:
        error_lines.append("=== STDOUT ===")
        error_lines.extend(result.stdout.split('\n')[-20:])  # Last 20 lines of stdout
    if result.stderr:
        error_lines.append("\n=== STDERR ===")
        error_lines.extend(result.stderr.split('\n')[-20:])  # Last 20 lines of stderr
    return "\n".join(error_lines)

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        os.chdir(os.path.join(working_dir, 'fastq', key, accession))

        success = False
        last_errors = {}

        for chemistry in chemistry_options:
            try:
                print(f"Trying chemistry {chemistry} for sample {accession}")
                cmd = [
                    "cellranger", "count", "--id="+accession+sample,
                    "--fastqs="+os.path.join(str(working_dir), 'fastq', str(key), str(accession)),
                    "--sample="+accession,
                    "--transcriptome="+os.path.join(os.environ['HOME'], 'WORKING/SC/ref/GRCh38'),
                    "--localcores="+str(NCPUS),
                    "--create-bam=true",
                    "--chemistry="+chemistry
                ]

                result = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )

                if result.returncode == 0:
                    print(f'Success with chemistry {chemistry} for sample {accession}')
                    successful_df = gse_dict[key][gse_dict[key]['run_accession'] == str(accession)]
                    successful_runs[key] = pd.concat([successful_runs[key], successful_df], ignore_index=True)
                    success = True
                    break
                else:
                    error_output = format_error_output(result)
                    last_errors[chemistry] = error_output
                    print(f'Chemistry {chemistry} failed for sample {accession}')
                    print(f"Error details:\n{error_output}\n")

            except Exception as e:
                error_msg = f"Exception with chemistry {chemistry}:\n{traceback.format_exc()}"
                last_errors[chemistry] = error_msg
                print(error_msg)
                continue

        if not success:
            print(f'All chemistry options failed for sample {accession}')
            error_msg = f"Exception with chemistry {chemistry}:\n{traceback.format_exc()}"
            last_errors[chemistry] = error_msg
            print(error_msg)

#%% Update gse_dict

metadata_dir = os.path.join(working_dir, 'metadata')
os.chdir(metadata_dir)


with open('dictionary_file_filtered.pkl', 'wb') as pkl_file:
    pickle.dump(successful_runs, pkl_file)
    print('Dictionary saved successfully.')

#%% Summary
for key in successful_runs:
    print(f"{key}: {len(successful_runs[key])} successful runs")

#%% End file
import sys
sys.exit()
