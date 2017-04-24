#!/usr/bin/env bash

printf "[ %s: kmer analysis with kmergenie ]\n" \
    "$(date)"

my_kmergenie=bin/kmergenie/kmergenie

# shellcheck disable=SC1091
source "src/sh/bash_header"
# shellcheck disable=SC1091
source "src/sh/io_parser"

# make outdir
outdir="$(dirname "${other_output}")"
if [[ ! -d "${outdir}" ]]; then
    mkdir -p "${outdir}"
fi

# set up kmergenie
read_file="${outdir}/read_file"
cat <<- _EOF_ > "${read_file}"
${input_fq}
_EOF_

# build command
cmd=( "${my_kmergenie}" "${read_file}" --diploid 
          -t "${max_cpus}" -o "${outdir}/histogram" )

shopt -s extglob
printf "Final command line: "
printf "%s " "${cmd[@]//+([[:blank:]])/ }"
printf "\n"
shopt -u extglob

# run kmergenie
log_file="${outdir}/kmergenie.log.txt"
srun --ntasks=1 --cpus-per-task="${max_cpus}" --exclusive \
    --output="${log_file}" "${cmd[@]}" &

printf "[ %s: Waiting for kmergenie to finish ]\n" "$(date)"
FAIL=0
fail_wait

# log metadata
metadata_file="${outdir}/METADATA.csv"

printf "[ %s: Logging metadata ]\n" "$(date)"
printf "metadata_file: %s\n" "${metadata_file}"
cat <<- _EOF_ > "${metadata_file}"
    Script,${0}
    branch,$(git rev-parse --abbrev-ref HEAD)
    hash,$(git rev-parse HEAD)
    date,$(date +%F)
    kmergenie version,$(my_kmergenie 2>&1 | sed '1q;d')
_EOF_

printf "[ %s: Done ]\n" "$(date)"

exit 0


