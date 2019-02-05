#!/bin/bash 
# POSIX

# Input parameters
data_file=''
max_parents=-1
max_tree_width=-1
twilp_output=/dev/null

while :; do
    case $1 in
        -f|--file)
            data_file=$2
            shift
            ;;
        -p|--parents)
            max_parents=$2
            shift
            ;;   
        -t|--tree_width)
            max_tree_width=$2    
            shift
            ;;
        -c|--clear)
            if [ -n "$(ls -A output)" ]; then
                $(rm output/*)
            fi
            if [ -n "$(ls -A logs)" ]; then
                $(rm logs/*)
            fi
            ;;
        -v|--verbose)
            twilp_output=/dev/stdout
            ;;
        --)
            shift
            break
            ;;
        *)
            break
    esac

    shift
done

if [ -z $data_file ]; then
    printf "No input file provided \n"; exit 1
fi

if [ ! -f $data_file ]; then
    printf "Input file doesn't exist \n"; exit 1
fi

if [ $max_parents -lt 0 -o $max_tree_width -lt 0 ]; then
    printf "Both max parents and max tree width must be set \n"; exit 1;
fi

# Strip the path and file extension from data_file
file_prefix=$(echo $data_file | cut -d '.' -f 1)
file_prefix=$(echo $file_prefix | awk -F'/' '{print $(NF)}')

# Running gobnilp score algorithm
score_file="output/$file_prefix.scores"
ESS=1

printf "Running gobnilp with ESS=$ESS and max parents=$max_parents \n"
./gobnilp1.2/scoring $data_file $ESS $max_parents > $score_file
printf "Score file stored to $score_file \n\n"

# Running TWILP
structure_output="output/"
max_time_pr_sub_ip=60
max_time_total_twilp=8000

printf "Running twilp with max tree width: $max_tree_width \n"
python ./twilp/twilp.py -f $score_file -o $structure_output -t $max_tree_width -p $max_parents -r $max_time_total_twilp -s $max_time_pr_sub_ip > "logs/twilp.log"
printf "BN structure stored to $structure_output \n\n"

# Running pgm.py 
## Parameter estimation
gml_file="${structure_output}tw_${max_tree_width}_mp_${max_parents}_${file_prefix}.scores_z.gml"
bif_output="output/$file_prefix.bif"

printf "Fitting the data.. \n"
source ./pgmpy/.myprojectenv/bin/activate 
python ./pgmpy/pgmpy-estimation.py -t $data_file -g $gml_file -o $bif_output > "logs/pgmpy-estimation.log"
printf "Saved to $bif_output \n\n"

## Inference algorithm
printf "Running inference.. \n"
# python ./pgmpy/pgmpy-inference.py -f $bif_output > "logs/pgmpy-inference.log"
printf "Done, results saved to output/$file_prefix.results \n\n"

## Calculate structural distances
printf "Calculating distances in true and learnt network \n"
network="alarm"
true_network="data/true_$network.bif"
translate_file="${network}_translate.txt"
python ./pgmpy/pgmpy-structural-distance.py -f "output/$file_prefix.bif" -t $true_network -tr $translate_file > "logs/pgmpy-structural-distances.log"
printf "Done, results of distance calculation saved to logs/pgmpy-structural-distances.log \n\n"

deactivate # deactivate virtual env for pgmpy