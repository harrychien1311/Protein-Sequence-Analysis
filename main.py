import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import os

def preprocess_tool_A(row):
    """Split protein names in one row of tool A result's dataframe 

    Args:
        row (_type_): a row of the dataframe

    Returns:
        dataframe: a preprocessed dataframe
    """
    proteins = row['ProteinName'].split(';')
    df = pd.DataFrame({'ProteinName': proteins, 'Sequence': [row['Sequence']] * len(proteins)})
    return df

def preprocess_tool_B(protein):
    """Extracting the real protein name of tool B result's dataframe

    Args:
        protein (str): protein name along with metadata information

    Returns:
        str: protein name with no metadata information
    """
    return protein.split("|")[2]

def extract_fasta(fasta_file):
    """Extract protein name and its sequence

    Args:
        fasta_file (str): the fastafile containing protein information
    Returns:
        dict: protein name and its sequence
    """
    protein_info = dict()
    # Parse the FASTA file
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Extract sequence ID and sequence
        protein_name = seq_record.id.split("|")[2]
        sequence = str(seq_record.seq)
        protein_info[protein_name] = sequence
    return protein_info

def is_simplified_semi_tryptic(row, protein_info):
    """Check if a sequene of peptide is a simplified semi tryptic

    Args:
        row: a row of result dataframe
        protein_info (dict): extracted protein and its sequence

    Returns:
        int: 1 if simplified semi tryptic otherwise 0
    """
    # Check if the protein name exists in the extracted file
    if not row["ProteinName"] in protein_info:
        return 0 
    # Extract sequence of the corresponding protein name
    sequence = protein_info[row["ProteinName"]]
    # Find the first index of the peptide in the sequence
    found_index = sequence.find(row["Sequence"])
    if found_index == -1:
        return 0
    elif found_index == 0 and (row["Sequence"].endswith("K") or row["Sequence"].endswith("R")):
        return 1
    elif found_index > 0:
        if sequence[found_index-1]!="K" and sequence[found_index-1]!="R" and (row["Sequence"].endswith("K") or row["Sequence"].endswith("R")):
            return 1
        elif sequence[found_index-1]=="K" or sequence[found_index-1]=="R" and not (row["Sequence"].endswith("K") or row["Sequence"].endswith("R")):
            return 1
        else: 
            return 0
    else:
        return 0
    
def calculate_protein_coverage(protein, tool_result, protein_info):
    """Calulate the total calculate protein coverage of each tool

    Args:
        protein (str): protein name
        tool_result (dataframe): result dataframe of tool A or tool B
        protein_info (dict): extracted protein information from fasta file

    Returns:
        float: the total protein coverage of tool's dataframe
    """
    peptides = list(tool_result[tool_result["ProteinName"]==protein]["Sequence"])
    if not protein in protein_info:
        return 0
    sequence = protein_info[protein]
    sequence_len = len(sequence)
    total_peptide_len = 0
    covered_indices = set()  # To keep track of the indices already covered by previous peptidess
    
    for peptide in peptides:
        index = sequence.find(peptide)
        
        # Check if the peptide is found and not overlapping with previously processed peptides
        if index != -1:
            overlapping = False
            for i in range(index, index + len(peptide)):
                if i in covered_indices:
                    overlapping = True
                    break
            if not overlapping:
                total_peptide_len += len(peptide)
                # Update covered indices
                covered_indices.update(range(index, index + len(peptide)))
            else:
                # If the peptide overlaps, calculate the length of non-overlapping part
                non_overlapping_length = len(peptide)
                for i in range(index, index + len(peptide)):
                    if i in covered_indices:
                        non_overlapping_length -= 1
                total_peptide_len += non_overlapping_length
        
    protein_coverage = total_peptide_len/sequence_len
    return protein_coverage
def generate_plot(values, titile, name, color):
    """Generate a specified bar plot showing the analyzed metrics
    for the report

    Args:
        values (list): value of calculated metrics
        titile (str): title of the plot
        name (str): name of the saved plot file
        color (str): color (specified in matplotlib) name for the bar plot
    """
    fig = plt.figure(figsize = (10, 5))
    plt.bar(["Tool A", "Tool B"], values, color=color, width=0.25)
    plt.title(titile, fontweight="bold")
    for i, value in enumerate(values):
        if name == "result_simplfied":
            plt.text(i, value+0.01, str(round(value,5)), ha='center', va='bottom')
            plt.ylim(0, 0.8)
        plt.text(i, value+0.5, str(round(value, 5)), ha='center', va='bottom')
    plt.grid(True)
    os.makedirs("figures", exist_ok=True)
    plt.savefig(f"figures/{name}.png")

if __name__=="__main__":
    # Read the data from tool A and tool B
    tool_A_result = pd.read_csv("tool_A_results.tsv", sep="\t")
    tool_B_result = pd.read_csv("tool_B_results.tsv", sep="\t")
    
    # Prepocess tool A
    tool_A_result.drop("index", axis=1, inplace=True)
    tool_A_result.rename(columns={"StrippedSequence": "Sequence"}, inplace=True)
    tool_A_result = pd.concat([preprocess_tool_A(row) for _, row in tool_A_result.iterrows()], ignore_index=True)
    # Preprocess tool B
    tool_B_result["ProteinName"] = tool_B_result["ProteinName"].apply(preprocess_tool_B)
     
    # Calculate number of different peptides
    nb_diff_peptide_tool_A = tool_A_result["Sequence"].nunique()
    nb_diff_peptide_tool_B = tool_B_result["Sequence"].nunique()
    
    # Calculate number of different proteins
    diff_protein_tool_A = list(tool_A_result["ProteinName"].unique())
    diff_protein_tool_B = list(tool_B_result["ProteinName"].unique())
    nb_diff_protein_tool_A = len(diff_protein_tool_A)
    nb_diff_protein_tool_B = len(diff_protein_tool_B)
    
    # Calculate number of Simplified Semi Tryptic
    protein_info = extract_fasta("uniprot_sprot_2021-01-04_HUMAN.fasta") # Extract the protein info from the fasta file
    tool_A_result["is_simplified"] = tool_A_result.apply(lambda row: is_simplified_semi_tryptic(row, protein_info), axis=1)
    tool_B_result["is_simplified"] = tool_B_result.apply(lambda row: is_simplified_semi_tryptic(row, protein_info), axis=1)
    nb_simpl_tool_A = len(tool_A_result[tool_A_result["is_simplified"]==1]) # Number of simplified semi tryptic in tool A
    ratio_simp_tool_A = nb_simpl_tool_A/len(tool_A_result["is_simplified"]) # Ratio of simplified semi tryptic in tool A
    nb_simpl_tool_B = len(tool_B_result[tool_B_result["is_simplified"]==1]) # Number of simplified semi tryptic in tool B
    ratio_simp_tool_B = nb_simpl_tool_B/len(tool_B_result["is_simplified"]) # Ratio of simplified semi tryptic in tool B
    
    #Calculate Protein coverage
    protein_coverage_tool_A = 0
    for protein in diff_protein_tool_A:
        protein_coverage_tool_A+= calculate_protein_coverage(protein, tool_A_result, protein_info)
    
    protein_coverage_tool_B = 0
    for protein in diff_protein_tool_B:
        protein_coverage_tool_B+= calculate_protein_coverage(protein, tool_B_result, protein_info)
    # Generate plots
    generate_plot([nb_diff_peptide_tool_A, nb_diff_peptide_tool_B], "Number of different peptides", "result_peptides", "r")
    generate_plot([nb_diff_protein_tool_A, nb_diff_protein_tool_B], "Number of different proteins", "result_proteins", "b")
    generate_plot([ratio_simp_tool_A, ratio_simp_tool_B], "Ratio of simplified-semi-tryptic peptides", "result_simplfied","g")
    generate_plot([protein_coverage_tool_A, protein_coverage_tool_B], "Protein sequence coverage", "result_protein_coverage", "orange")
    