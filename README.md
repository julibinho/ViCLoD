# ViCloD

**VIsualizing BCR repertoires and intraCLonal Diversity**

ViCloD is a web server for visualizingz the clonality and intraclonal diversity of BCR repertoires. Here we make the ViCloD pipeline available, the stand-alone version. Users can run it locally and visualize their results online. This pipeline uses the following tools: [MobiLLe](https://github.com/julibinho/MobiLLe) to cluster the clonally-related sequences in BCR repertoires and [ClonalTree](https://github.com/julibinho/ClonalTree) to reconstructing the evolutionary history of a BCR lineages.

**REFERENCE**
Lucile Jeusset, Nika Abdollahi, Anne Langlois De Septenville, Marine Armand, Thibaud Verny, Clotilde Bavetti, Frédéric Davi and Juliana S. Bernardes. ViCloD, an interactive web tool for visualizing B cell repertoires and analyzing intra-clonal diversities in B-cell tumors. To be submitted.

**CONTACT**
  E-mail: 
  juliana.silva_bernardes@sorbonne-universite.fr 
  
## Inputs

  * An AIRR formatted file containing annotated IGH sequences.
  * See [example input file](https://github.com/julibinho/ViCLoD/Examples/Input)

## Outputs

  * ViCloD returns:
  
    - [repertoire_name]_visualization.zip : contains all the files necessary for viewing the BCR repertoire in the format requested by the ViCloD website (compressed file)

    - [repertoire_name]_repertoire.json : the clustering output of clonally-related sequences of the repertoire in [JSON format](https://fr.wikipedia.org/wiki/JavaScript_Object_Notation)
    
    - [repertoire_name]_tree_all.json : the reconstructed BCR lineage tree of the 5 first clones in [JSON format](https://fr.wikipedia.org/wiki/JavaScript_Object_Notation) 
    
    - [repertoire_name]_tree_simplification1.json : the first level of simplification of the previously reconstructed BCR lineage tree in [JSON format](https://fr.wikipedia.org/wiki/JavaScript_Object_Notation)
    
    - [repertoire_name]_tree_simplification2.json : the second level of simplification of the previously reconstructed BCR lineage tree in [JSON format](https://fr.wikipedia.org/wiki/JavaScript_Object_Notation)

    - [repertoire_name]_C[clone_number]_HAUS_sequence.fasta : a fasta file containing the HAUS sequence of the indicated clone
    
    - [repertoire_name]_C[clone_number]_subclones_sequences.fasta : a fasta file containing the subclone sequences of the indicated clone
    
    - [repertoire_name]_log.txt : contain the information on the reads used for the analysis.
    
    - [repertoire_name]_unannotated_seq.txt : contains all the sequences of the input AIRR file which are not correctly annotated to allow analysis by ViCloD
    
    - [repertoire_name]_AIRR_seq_name_matching.txt : contain the match between the original name of the reads (in AIRR file) and the new id provide by ViCloD
    
    - [repertoire_name]_unique_seq_id_matching.txt : ids of sequences with the same nucleotide sequence
  * See [example output files](https://github.com/julibinho/ViCLoD/Examples/Output)
     
      
## Requirements 

  * We strongly recommend [anaconda](https://docs.anaconda.com/anaconda/install/) environment. 
  
  * Python version 3 or later

  * numpy :
      ```
      conda install numpy
      ```
      or 
      ```
      pip install numpy
      ```

  * matplotlib
    ```
      conda install -c conda-forge matplotlib
     ```
     or
  
      ```
      pip install matplotlib
      ```
      
  * Palettable :
      ```
      conda install -c conda-forge palettable
      ```
      or
      ```
      pip install palettable
      ```

  * skbio
      ```
      conda install -c anaconda scikit-bio
      ```
      or 
      ```
      pip install scikit-bio
      ```
  * Levenshtein
      ```
      conda install -c conda-forge python-levenshtein 
      ```
      or
      ```
      pip install python-Levenshtein
      ```

  * Biopython
      ```
      conda install -c conda-forge biopython
      ```
      or
      ```
      pip install biopython
      ```

  * ete3 :
      ```
      conda install -c etetoolkit ete3
      ```
      or
      ```
      pip install ete3
      ```
  
  * python-newick :
      ```
      conda install -c bioconda python-newick
      ```
      or
      ```
      pip install newick
      ```
  
  * muscle :
      ```
      conda install -c bioconda muscle
      ```
      or
      ```
      pip install muscle
      ```

## Using ViCloD 
   The command line for launching the ViCloD is:

  ```
  $ bash run_ViCloD.sh [AIRR_file] [output_path] [threshold_value] [thresold_type]

  ```
### required arguments 
  * [AIRR_file] is the AIRR file containing the annotated sequences of the BCR repertoire,
  * [output_path] is path of the ouptut directory
  * [threshold_value] is a threshold for eliminating infrequent reads
  * [thresold_type] is the type in which the threshold value is given (percentage or number)


  For instance the following command can be run in the src/ folder:
  ```
  $ src/run_ViCloD.sh Examples/Input/example.tsv Examples/Output 0 number
  ```
                      
  Output files will be placed in a folder as such:
  ```
  [output_path]/[AIRR_file_name]_ViCloD_output
  ```
  * [AIRR_file_name] is the name of the AIRR file provide as input
  
## License, Patches, and Ongoing Developements

  * The program is distributed under the CeCILL licence 
  * [Feature requests and open issues](https://github.com/julibinho/ViCloD/issues).
 
 
