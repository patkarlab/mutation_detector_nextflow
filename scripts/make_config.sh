PWD="/home/pipelines/mutation_detector_nextflow"
while getopts s: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
    esac
done

BAM=$PWD/$sample"/gatk38_processing"/$sample".final.bam"
echo $sample
echo $PWD
echo $BAM
printf '%s\t%s\t%s' "$BAM" "600" "$sample" > $PWD/config.txt
