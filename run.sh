#!/bin/bash

#######################################
### Run spatial clustering analysis ###
#######################################

function usage () {
    cat >&2 <<EOF

USAGE: ${0} [-y] [options]
  -y <config file> : Path to the YAML config file. Required.
  -s Submit job.
  -v Verbose.
  -h Print the usage info.

EOF
}
# initial : makes this loop silent and now requires '?)'
# ${opt} is each option and ${OPTARG} its the argumet (if a colon is there ${opt}:)
SUBMIT=FALSE
VERBOSE=TRUE
while getopts ":y:s:vh" opt; do
  case ${opt} in
    y) CONFIG_FILE=${OPTARG};;
    s) SUBMIT=TRUE;;
    v) VERBOSE=TRUE;;
    h) usage; exit 1;;
    \?) echo "No -${OPTARG} argument found."; usage; exit 1;;
  esac
done
if [[ ${OPTIND} -eq 1 ]] ; then
    usage; exit 1
fi

function read_yaml(){
  sed 's/#.*//g' ${1} | grep ${2}: | sed 's/.*:[^:\/\/]//; s/\"//g'
}
# CONFIG_FILE="/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/spatial_pilot_5k/scripts/devel/config.yaml"
# CONFIG_FILE="/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/spatial_pilot_5k/scripts/devel/PBT_config.yaml"

# CONFIG_FILE=/home/fcastaneda/bin/spatial_clustering/config.yaml

OUTPUT_DIR="$(read_yaml ${CONFIG_FILE} output_dir)"
PROJECT_NAME="$(read_yaml ${CONFIG_FILE} project_name)"
OUTPUT_DIR="${OUTPUT_DIR%/}/${PROJECT_NAME}"
PIPELINE_FOLDER="$(read_yaml ${CONFIG_FILE} pipeline)"
CONDA_ENV="$(read_yaml ${CONFIG_FILE} conda_env)"
R_MODULE="$(read_yaml ${CONFIG_FILE} R_module)"
CLUSTER_CONFIG="$(read_yaml ${CONFIG_FILE} cluster_config)"
PROJECT_NAME="$(read_yaml ${CONFIG_FILE} project_name)"


if [[ -v R_MODULE ]]; then module load $R_MODULE; else echo "Using default R"; fi

echo -e "\033[0;36m------------------------------- PRESENTING PARAMETERS -------------------------------\033[0m"
echo "Configuration file: ${CONFIG_FILE}"
echo "Output path: ${OUTPUT_DIR}"

if [[ ! -d "${OUTPUT_DIR}" ]]; then mkdir --parents "${OUTPUT_DIR}"; fi
if [[ ! -d "${OUTPUT_DIR}/scripts" ]]; then mkdir "${OUTPUT_DIR}/scripts"; fi
cd ${OUTPUT_DIR}

# rm -r "${OUTPUT_DIR}/scripts/*" # Is this necessary?

cp -r ${PIPELINE_FOLDER}/slurm .

if [[ $(wc -l $CLUSTER_CONFIG) > 0 ]]; then 
  sed -i 's|/home/fcastaneda/bin/spatial_clustering/cluster_slurm.json|'"${CLUSTER_CONFIG}"'|g' ${OUTPUT_DIR}/slurm/config.yaml;
  echo "Using user defined HPC resources";
  else echo "Using HPC resouces defined in the slurm folder"; 
fi

cp $CONFIG_FILE config.yaml
# conda activate $CONDA_ENV

JOBFILE="${OUTPUT_DIR}/scripts/clump_${PROJECT_NAME}"

cp /home/fcastaneda/bin/spatial_clustering/routine_template_emma_Slurm_Spatial_Clustering.sh ${JOBFILE}.sh # for the new HPC

sed -i 's|{cellranger}|spatial|' ${JOBFILE}.sh
sed -i 's|{username}|'"${USER}"'|g' ${JOBFILE}.sh
sed -i 's|{sampleid}|'"${PROJECT_NAME}"'|g' ${JOBFILE}.sh
sed -i 's|\/\.\.||g' ${JOBFILE}.sh
sed -i 's|beegfs|BioScratch|g' ${JOBFILE}.sh #Change due problem in beegfs
sed -i 's|{routine_pbs}|clump|' ${JOBFILE}.sh
sed -i 's|{outpath}|'"${OUTPUT_DIR}"'|g' ${JOBFILE}.sh
sed -i 's|{PIPELINE_FOLDER}|'"${PIPELINE_FOLDER}"'|g' ${JOBFILE}.sh
# sed -i 's|{OUTPUT_DIR}|'"${OUTPUT_DIR}"'|g' ${JOBFILE}.sh

if [[ -v R_MODULE ]]; then sed -i 's|#module load {R_module}|module load '"${R_MODULE}"'|g' ${JOBFILE}.sh; echo "Using user defined R"; else echo "Using default R"; fi

sed -i 's|#conda activate {conda_env}|conda activate '"${CONDA_ENV}"'|g' ${JOBFILE}.sh

sed -i 's|cp ${PROJ.*|cp -r ${PROJDIR}/. ./|g' ${JOBFILE}.sh # to copy everything to scratch
sed -i 's|cp -R ./.*${PROJ.*|cp -r . ${PROJDIR}/|g' ${JOBFILE}.sh # copy from scratch

echo "Job file: ${JOBFILE}.sh"
echo "Pushing critical lines..."

sed -i 's|{walltime}|100:00:00|g' ${JOBFILE}.sh
sed -i 's|{nodes}|1|g' ${JOBFILE}.sh
sed -i 's|{ppn}|1|g' ${JOBFILE}.sh
sed -i 's|{mem}|8gb|g' ${JOBFILE}.sh

if echo "${SUBMIT}" | grep -qE "TRUE|^yes$|^y$"; then
  echo "Check it out"; exit
fi
sbatch ${JOBFILE}.sh


 