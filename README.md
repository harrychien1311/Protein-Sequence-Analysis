# Protein Sequence extraction result analysis
In this repo, I analyzed the protein sequence extracted by using 2 tools Javelin and Skeleton-productions of [Bionsight](https://www.bionsight.com/platform/).
 Based on the analyzed result I can conclude which tool gave a better performance in finding protein sequence.
 The detail of the requirements is viewed [here](https://bionsight.notion.site/Basics-of-Proteomics-DS-69a01d4779404af6bf10f06d673fcfb8)

## Criteria

The measures for “better performance” are as follows:

1. Which tool found more different `Peptides`
2. Which tool found more different `Proteins`
3. When the ratio of `Simplified-Semi-Tryptic` among the found `Peptides` is high
     1. The case of `Simplified-Semi-Tryptic` is as follows.
         1. If the `Sequence` of `Peptide` ends with K or R, and the amino acid immediately before the `Sequence` of `Peptide` in the `Sequence` of `Protein` does not end with K or R.
         2. If the `Sequence` of `Peptide` does not end with K or R, and the amino acid immediately before the `Sequence` of `Peptide` in the `Sequence` of `Protein` ends with K or R.
        
         <aside>
         💡 1. If the `Sequence` of `Peptide` begins at the beginning of the `Sequence` of `Protein` (if there is no preceding amino acid), it is considered to end with K or R.
         2. If `Sequence` of `Peptide` appears multiple times in `Sequence` of `Protein`, it shall be the first occurrence.
         3. If there is no `Protein` name information, it is not `Simplified-Semi-Tryptic`
        
         </aside>
        
     2. An example of `Simplified-Semi-Tryptic` is [below](https://www.notion.so/Basics-of-Proteomics-DS-69a01d4779404af6bf10f06d673fcfb8?pvs=21)
4. When `Protein sequence coverage` is high
    
     <aside>
     💡 Protein sequence coverage of `Protein A` is as follows
     [Example below](https://www.notion.so/Basics-of-Proteomics-DS-69a01d4779404af6bf10f06d673fcfb8?pvs=21)
     </aside>

## Installation
Create new environment with python=3.11
```bash
conda create -n assigment python=3.11
```

Install dependencies:
```bash
pip install -r requirements.txt
```
### How to run
```bash
python3 main.py
```