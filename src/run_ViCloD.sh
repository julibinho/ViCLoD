#!/bin/bash

BASEDIR=$(dirname "${0}")
cd "$BASEDIR"

#check the number of argument
if [[ $# -ne 2 ]]; then
  echo "Invalid number of arguments"
  exit
fi

filename=$1
extension="${filename##*.}"
path=$(dirname -- "$filename")
analysis=$2

#check the format of the file
if [ "$extension" != "tsv" ] && [ "$extension" != "txt" ]; then
  echo "Please provide a file in AIRR format."
  exit
fi

#check if all the information are in the AIRR file
header=$(head -n 1 "$filename" | grep -w -iEo "sequence_id|sequence|productive|d_call|j_call|fwr1|cdr1|fwr2|cdr2|fwr3|cdr3|fwr4|v_germline_alignment|np1|d_germline_alignment|np2|j_germline_alignment|cdr1_start|cdr2_start|cdr3_start|cdr1_end|cdr2_end|cdr3_end|v_sequence_start" | sort | uniq -c | wc -l)
if [[ $header -ne 24 ]]; then
  echo "Information are missing in the AIRR file, please provide a file in the correct format."
  exit
fi

#check if the output file has already been created
if [[ -d "${path}/${analysis}_ViCloD_output" ]]
then
  echo "${analysis}_ViCloD_output exists, please remove it before continuing."
  exit
fi

mkdir ${path}/${analysis}_ViCloD_output
output=${path}/${analysis}_ViCloD_output

BEFORE=$SECONDS

echo "## MobiLLe analysis :"

python MobiLLe/check_airr_file_informations.py -a ${filename} -o ${output}/${analysis}
python MobiLLe/format_labeling_imgt_airr.py -a ${output}/${analysis}_new_id.tsv -o ${output}/${analysis} -t 0.005 
python MobiLLe/initial_clustering_f.py -i ${output}/${analysis}_seq_Fo_V_CDR3_Jseq.txt -o ${output}/${analysis} -s 0.7 
python MobiLLe/format_clustering_output.py -i ${output}/${analysis}_sameVJ_noallele_CDR3_0.7.txt -o ${output}/${analysis}_initial_clusters_Fo.txt 
python MobiLLe/refinement.py -f ${output}/${analysis}_seq_Fo_V_CDR3_Jseq.txt -c ${output}/${analysis}_initial_clusters_Fo.txt
python MobiLLe/format_clustering_output.py -i ${output}/${analysis}_seq_Fo_V_CDR3_Jseq_clone_V_CDR3_J.txt -o ${output}/${analysis}_final_clusters_Fo.txt 
python MobiLLe/two_level_clonal_info.py -f ${output}/${analysis}_seq_Fo_V_CDR3_Jseq.txt -c ${output}/${analysis}_final_clusters_Fo.txt -n ${output}/${analysis} 

python select_clone/intraclonal_study_input_creator.py -r ${output}/${analysis}_repertoire_two_levels_info.txt -c ${output}/${analysis}_final_clusters_seq_info.txt -n ${output}/${analysis} -s 5 

echo "## Tree generation for the 5 most abundant clonotypes with clonalTree :"
for i in {1..5}
  do
    if [ -e ${output}/${analysis}_${i}_seq_info.txt ]
    then
      echo "  # Clone ${i}"
      python select_clone/alignment_intraclonal.py -a ${output}/${analysis}_new_id.tsv -f ${output}/${analysis}_${i}_seq_info.txt -n ${output}/${analysis}_${i} -s 200 
      python clonalTree/clonalTree.py  -i ${output}/${analysis}_${i}_200_selected_seq_uniq.aln.fa -o ${output}/${analysis}_${i}_200.nk -a 0 -r 0 
      python clonalTree/SimplificationVF.py -n ${output}/${analysis}_${i}_200.nk -a ${output}/${analysis}_${i}_200_selected_seq_uniq.aln.fa
    else
      echo "  # There is no clone C${i}"
    fi
  done
  
echo "## Visualization :"
python visualization/clone_informations_json.py -r ${output}/${analysis}_repertoire_two_levels_info.txt -l ${output}/${analysis}_log.txt -a ${analysis} -c visualization/colors.txt
python visualization/generate_clonotype_files.py -c 5 -o ${output}/ -a ${analysis} -n _200.nk -t all
python visualization/generate_clonotype_files.py -c 5 -o ${output}/ -a ${analysis} -n _200_simplification1.nk -t simplification1
python visualization/generate_clonotype_files.py -c 5 -o ${output}/ -a ${analysis} -n _200_simplification2.nk -t simplification2

#rm ${output}/${analysis}_*_clusters_* ${output}/${analysis}_sameVJ_noallele_CDR3_0.7.txt ${output}/${analysis}_seq_Fo_V_CDR3_* ${output}/${analysis}_repertoire_two_levels_info.txt
#rm ${output}/${analysis}*_200* ${output}/${analysis}*_seq_info.txt

zip -j ${output}/${analysis}_visualization.zip ${output}/${analysis}_*.json ${output}/${analysis}_C[0-9]_simplification2_sequences.txt

ELAPSED=$(($SECONDS-$BEFORE))
echo "ViCloD pipeline run time : ${ELAPSED} seconds"

