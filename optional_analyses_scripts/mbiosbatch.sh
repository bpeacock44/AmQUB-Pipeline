#!/bin/bash
# Usage: mbiosbatch -S <step> -J <jobname> [-T <time>] [-C <cpus>] [-M <mem>] [-P <partition>] [-- pipeline args]
# Example: mbiosbatch -S 1 -J JB141_part1 -T 04:00:00 -C 8 -M 16G -P short -- -f data/JB141.fq -p data/JB141_map.txt -m 3
STEP=""
JOBNAME="mbio_job"
TIME=""
CPUS="32"
MEM=""
PARTITION=""
EMAIL="beth.b.peacock@gmail.com"

while getopts ":S:J:T:C:M:P:E:" opt; do
    case $opt in
        S) STEP="$OPTARG" ;;
        J) JOBNAME="$OPTARG" ;;
        T) TIME="$OPTARG" ;;
        C) CPUS="$OPTARG" ;;
        M) MEM="$OPTARG" ;;
        P) PARTITION="$OPTARG" ;;
        E) EMAIL="$OPTARG" ;;
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
${PARTITION:+#SBATCH --partition=${PARTITION}}
#SBATCH --cpus-per-task=${CPUS}
${MEM:+#SBATCH --mem=${MEM}}
${TIME:+#SBATCH --time=${TIME}}
#SBATCH --output=${JOBNAME}_%j.out
#SBATCH --error=${JOBNAME}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}

module load db-ncbi
module load singularity
singularity exec \\
    --bind ${programs}:/home/programs/ \\
    --bind \${NCBI_DB}:/database/ \\
    --bind \$(pwd):/home/analysis/ \\
    ${programs}/testing.sif \\
    AmQUB_part${STEP}.sh $@
EOF