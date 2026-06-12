#!/bin/bash
# Usage: mbiosbatch -s <step> -j <jobname> [-t <time>] [-c <cpus>] [-m <mem>] [-p <partition>] [-- pipeline args]
# Example: mbiosbatch -s 1 -j JB141_part1 -t 04:00:00 -c 8 -m 16G -p short -- -f data/JB141.fq -p data/JB141_map.txt -m 3

STEP=""
JOBNAME="mbio_job"
TIME="7-00:00:00"
CPUS="32"
MEM="32G"
PARTITION="epyc"

while getopts ":s:j:t:c:m:p:" opt; do
    case $opt in
        s) STEP="$OPTARG" ;;
        j) JOBNAME="$OPTARG" ;;
        t) TIME="$OPTARG" ;;
        c) CPUS="$OPTARG" ;;
        m) MEM="$OPTARG" ;;
        p) PARTITION="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
[[ "$1" == "--" ]] && shift

if [[ -z "$STEP" ]]; then
    echo "Error: -s <step> is required." >&2
    exit 1
fi

programs="/rhome/bpeacock/shared/mbio_pipeline_files"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${JOBNAME}
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --output=${JOBNAME}_%j.out
#SBATCH --error=${JOBNAME}_%j.err

module load db-ncbi
module load singularity

singularity exec \\
    --bind ${programs}:/home/programs/ \\
    --bind \${NCBI_DB}:/database/ \\
    --bind \$(pwd):/home/analysis/ \\
    ${programs}/testing.sif \\
    /home/programs/AmQUB_part${STEP}.sh $@
EOF