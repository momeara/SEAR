Lazarus of Bethany was brought back to life four days after his death, as described int he Gospel of John.

Lazarus is also the name of the Thornton lab program for ancestral
gene reconstruction. This R package provides context and tools to help
use Lazarus. Including, automatically building trees form sets of
sequences.

## to install:
1. install required software and make sure it's on your path
   * blast+
   * clustalo
   * readseq
   * phyml

2) From within R:

    install.packages("devtools")
    devtools::install_github("momeara/Bethany")

# usage
1. Use http://www.ebi.ac.uk/Tools/hmmer to retrieve a fasta file. For PhyML, aim to get ~ 100 sequences or less.
2. Run pipeline
```    
sequences_to_trees(
    sequences = <path/to/sequences.fa.gz>,
    output_path = <prot_dir/data>,
    run_id = <query_tool_seqdb>,
    user_agent("httr <name@example.com>")) # for retrieving info from uniprot
```

This will produce a directory
```
prod_dir/data/<run_id>_YYMMDD/
   <sequences.fa.gz>
   # input fasta sequences

   run.sh
   # bash commands to re-run analysis
 
   name_map.xlsx
   # mapping between input, phylip and tree names, along with
   # the sequences an and information about each sequence
   # retrieved form Uniprot

   clustalo_input.fasta
   # input sequecences in fasta format with phylip names for clustal Omega
 
   clustalo_output.phylip3.2
   # multiple sequence alignment with phylip names in phylip3.2 format
 
   phyml_input.phylip
   # multiple sequence alignment  in phylip4 format (phylip names)
 
   phyml_input.phylip_phyml_stats_<run_id>.txt
   # output log of running phyml
 
   phyml_input.phylip_phyml_tree_<run_id>.txt
   # tree generate as output from phyml (phylip names)
 
   <run_id>.newick
   # tree with informative names
```