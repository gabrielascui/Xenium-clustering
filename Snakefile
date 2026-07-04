"""
### Snakefile for Spatial transcriptomics
"""
import pandas as pd
import os
import re

conf_file_loc = "config.yaml"



configfile:conf_file_loc

# file_path = config["input_expression"] +  "transcripts.csv.gz"

# if os.path.exists(file_path):
#     print(f"Transcripts file '{file_path}' exists.")
# else:
#     print(f"Transcripts file don't exist. This will take som extra time '{file_path}'")

get_transcripts_csv_script=config["pipeline"] + "/R/get_transcripts_csv.R"
getting_expr_data=config["pipeline"] + "/R/getting_expression_data.R"
init_script=config["pipeline"] + "/R/init.R"
markers_script=config["pipeline"] + "/R/markers.R"
compo_script=config["pipeline"] + "/R/report_components.R"

MEANS = config['variable_features']['mean.cutoff'][0]
PERCENTAGES = config['variable_features']['percent']
COMPONENTES = config['dim_reduction']['base']['chosen_comp']
RESOLUTIONS = config['resolution']
EXECS = config['exec']

inputs_expression= config["input_expression"]
if isinstance(inputs_expression, str):
    inputs_expression = [inputs_expression]
inputs_expression_files = [i for i in range(1, len(inputs_expression)+1)]

rule all:
    input:
        expand("transcripts_done_slide{slide}.txt", slide = inputs_expression_files),
        ".object_input_data.rds",
        ".object_normalized_base.rds",
        expand(".object_norm_mean{mean}_pct{percentage}.rds", mean = MEANS, percentage = PERCENTAGES),
        expand(".object_pca_mean{mean}_pct{percentage}.rds", mean = MEANS, percentage = PERCENTAGES),
        expand(".object_init_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
        expand(".object_metadata_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
        expand(".object_reductions_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
        expand(".object_graphs_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
        expand(".markers_mean{mean}_pct{percentage}_pc{component}.txt", mean = MEANS ,percentage = PERCENTAGES, component = COMPONENTES),
        expand(".report_mean{mean}_pct{percentage}_pc{component}.txt", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES)

 ### --------------- Cellranger parquet file prepocessing --------------- ###
rule get_transcripts_csv:
    input:
        conf_file_loc
    output:
        "transcripts_done_slide{slide}.txt"
    shell: 
        '{EXECS} {get_transcripts_csv_script} --yaml {conf_file_loc} --slide {wildcards.slide} -v TRUE'

 ### --------------- Seurat Normalization and HVG selection --------------- ###

rule getting_expression_data:
    input:
        expand("transcripts_done_slide{slide}.txt", slide = inputs_expression_files)
    output:
        ".object_input_data.rds",    
    params:
        component = config['dim_reduction']['base']['n_comp']
    shell:
        "{EXECS} {getting_expr_data} -y {conf_file_loc}"


rule normalize_base:
    input:
        ".object_input_data.rds"
    output:
        ".object_normalized_base.rds"
    shell:
        "{EXECS} {init_script} -y {conf_file_loc} --stage normalize_base --prefix normalized_base"

rule hvg_object:
    input:
        ".object_normalized_base.rds"
    output:
        ".object_norm_mean{mean}_pct{percentage}.rds"
    shell:
        "{EXECS} {init_script} -y {conf_file_loc} --stage hvg --percent {wildcards.percentage} --prefix norm_mean{wildcards.mean}_pct{wildcards.percentage}"

rule pca_object:
    input:
        ".object_norm_mean{mean}_pct{percentage}.rds"
    output:
        ".object_pca_mean{mean}_pct{percentage}.rds"
    params:
        component = config['dim_reduction']['base']['n_comp']
    shell:
        "{EXECS} {init_script} -y {conf_file_loc} --stage pca --percent {wildcards.percentage} --n_comp {params.component} --prefix pca_mean{wildcards.mean}_pct{wildcards.percentage}"

rule init_object:
    input:
        ".object_pca_mean{mean}_pct{percentage}.rds"
    output:
        ".object_init_mean{mean}_pct{percentage}_pc{component}.rds",
        ".object_metadata_mean{mean}_pct{percentage}_pc{component}.rds",
        ".object_reductions_mean{mean}_pct{percentage}_pc{component}.rds",
        ".object_graphs_mean{mean}_pct{percentage}_pc{component}.rds"
    shell:
        "{EXECS} {init_script} -y {conf_file_loc} --stage cluster --percent {wildcards.percentage} --chosen_comp {wildcards.component} --prefix init_mean{wildcards.mean}_pct{wildcards.percentage}_pc{wildcards.component}"

rule markers:
    input:
        ".object_init_mean{mean}_pct{percentage}_pc{component}.rds"
    output:
        ".markers_mean{mean}_pct{percentage}_pc{component}.txt"
    message: " --- Branch resolution for marker calculation --- "
    shell:
        "{EXECS} {markers_script} -y {conf_file_loc} --init_file {input} --percent {wildcards.percentage} --chosen_comp {wildcards.component} --prefix init_mean{wildcards.mean}_pct{wildcards.percentage}_pc{wildcards.component}"

rule report_components:
    input:
        ".object_init_mean{mean}_pct{percentage}_pc{component}.rds"
    output:
        ".report_mean{mean}_pct{percentage}_pc{component}.txt"
    message: " --- Creating report: components  ---"
    shell:
        "{EXECS} {compo_script} -y {conf_file_loc} --init_file {input} --percent {wildcards.percentage} --chosen_comp {wildcards.component} --prefix init_mean{wildcards.mean}_pct{wildcards.percentage}_pc{wildcards.component}"
