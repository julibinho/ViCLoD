#!/bin/bash

BASEDIR=$(dirname "${0}")

#check the number of argument
if [[ $# -ne 4 ]]; then
  echo "Invalid number of arguments"
  exit
fi

input=$1
filename=$(basename -- "$input")
extension="${filename##*.}"
analysis="${filename%.*}"
path=$2
threshold=$3
threshold_type=$4
nb_clone=5
nb_sub_clone=200
gene_allele=0

#check if the input file exist
if [[ ! -f "${input}" ]]; then
  echo "The file passed as input doesn't exist. Please provide a correct file."
  exit
fi

#check the format of the file
if [ "$extension" != "tsv" ] && [ "$extension" != "txt" ]; then
  echo "Please provide a file in AIRR format."
  exit
fi

#check if the output file has already been created
if [[ -d "${path}/${analysis}_ViCloD_output" ]]
then
  echo "${analysis}_ViCloD_output exists, please remove it before continuing."
  exit
fi

mkdir -p ${path}/${analysis}_ViCloD_output
output=${path}/${analysis}_ViCloD_output

BEFORE=$SECONDS

echo "## Merge same lines in AIRR file"
header=$(head -n 1 ${input})
id_pos=$(echo -e "$header"  | awk -v RS='\t' '/sequence_id/ { print NR; exit }')
echo -e "sequence_old_ids	$header" > ${output}/${analysis}_vquest_airr.tsv
awk -v "elem=$id_pos" 'BEGIN{FS=OFS="\t"} FNR>1 { val=$elem;$elem="";a[$0]=a[$0]val";"; }END{ for (i in a)print a[i], i; }' ${input} >> ${output}/${analysis}_vquest_airr.tsv

echo "## Filter data from AIRR file"
python ${BASEDIR}/MobiLLe/format_airr_file_data.py -a ${output}/${analysis}_vquest_airr.tsv -o ${output}/${analysis} -m $threshold -t $threshold_type -q True
rm ${output}/${analysis}_vquest_airr.tsv

if [ -f "${output}/${analysis}_airr_new_id.tsv" ]; then

  echo "## MobiLLe analysis"
  python ${BASEDIR}/MobiLLe/format_labeling_imgt_airr.py -a ${output}/${analysis}_airr_new_id.tsv -o ${output}/${analysis} -g $gene_allele
  python ${BASEDIR}/MobiLLe/initial_clustering_f.py -i ${output}/${analysis}_seq_Fo_V_CDR3_Jseq.txt -o ${output}/${analysis} -s 0.7 -a ${output}/
  python ${BASEDIR}/MobiLLe/format_clustering_output.py -i ${output}/${analysis}_sameVJ_noallele_CDR3_0.7.txt -o ${output}/${analysis}_initial_clusters_Fo.txt
  python ${BASEDIR}/MobiLLe/refinement.py -i ${output}/${analysis}_seq_Fo_V_CDR3_Jseq.txt -l ${output}/${analysis}_initial_clusters_Fo.txt -v 1 -j 2 -c 2 -d 1 -m 0 -o ${output}/${analysis}
  python ${BASEDIR}/MobiLLe/format_clustering_output.py -i ${output}/${analysis}_clone_V_CDR3_J.txt -o ${output}/${analysis}_final_clusters_Fo.txt
  python ${BASEDIR}/MobiLLe/two_levels_info.py -a ${output}/${analysis}_airr_new_id.tsv -c ${output}/${analysis}_final_clusters_Fo.txt -o ${output}/${analysis} -g $gene_allele

  python ${BASEDIR}/select_clone/alignment_intraclonal.py -r ${output}/${analysis}_repertoire_two_levels_info.txt -c $nb_clone -s $nb_sub_clone -o ${output}/${analysis}

  echo "## Tree generation for the 5 most abundant subclones with clonalTree :"
  for (( i=1; i <= ${nb_clone}; i++ ))
  do
    if [ -e ${output}/${analysis}_${i}_${nb_sub_clone}_selected_seq_uniq.aln.fa ]
    then
      echo "  Clone ${i}"
      python ${BASEDIR}/clonalTree/clonalTree.py -i ${output}/${analysis}_${i}_${nb_sub_clone}_selected_seq_uniq.aln.fa -o ${output}/${analysis}_${i}_${nb_sub_clone}.nk -a 0 -r 0
      python ${BASEDIR}/clonalTree/SimplificationVF.py -n ${output}/${analysis}_${i}_${nb_sub_clone}.nk -a ${output}/${analysis}_${i}_${nb_sub_clone}_selected_seq_uniq.aln.fa
    else
      echo "  There is no clone C${i}"
    fi
  done

  echo "## Visualization :"
  python ${BASEDIR}/visualization/clone_information_json.py -r ${output}/${analysis}_repertoire_two_levels_info.txt -l ${output}/${analysis}_log.txt -p ${output}/${analysis}_clonotypes_proportion.txt -a ${analysis} -c ${BASEDIR}/visualization/colors.txt
  python ${BASEDIR}/visualization/subclone_information_json.py -c $nb_clone -s $nb_sub_clone -o ${output} -a ${analysis} -n ${nb_sub_clone}.nk -t all
  python ${BASEDIR}/visualization/subclone_information_json.py -c $nb_clone -s $nb_sub_clone -o ${output} -a ${analysis} -n ${nb_sub_clone}_simplification1.nk -t simplification1
  python ${BASEDIR}/visualization/subclone_information_json.py -c $nb_clone -s $nb_sub_clone -o ${output} -a ${analysis} -n ${nb_sub_clone}_simplification2.nk -t simplification2

  rm ${output}/${analysis}_*_clusters_* ${output}/${analysis}_sameVJ_noallele_CDR3_0.7.txt ${output}/${analysis}_seq_Fo_V_CDR3_* ${output}/${analysis}_clone_V_CDR3_J.txt ${output}/${analysis}_repertoire_two_levels_info.txt ${output}/${analysis}*_airr_new_id.tsv
  rm ${output}/${analysis}_clonotypes_proportion.txt ${output}/${analysis}*_200*
  
  zip -j ${output}/${analysis}_visualization.zip ${output}/${analysis}_*.json ${output}/${analysis}_*_HAUS_sequence.fasta ${output}/${analysis}_*_subclones_sequences.fasta ${output}/${analysis}_log.txt

fi

ELAPSED=$(($SECONDS-$BEFORE))
echo "Pipeline run time : ${ELAPSED} seconds"


