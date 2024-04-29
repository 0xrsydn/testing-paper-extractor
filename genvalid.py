from Bio import Entrez, Seq
import pandas as pd 

# Set up Entrez
Entrez.email = "rafayrussani@gmail.com"  # Provide your email for NCBI

def is_valid_gene(gene_name):
    try:
        # Query NCBI for the gene name
        handle = Entrez.esearch(db="gene", term=gene_name)
        record = Entrez.read(handle)
        if record['IdList']:
            return True
        else:
            # Try fuzzy matching for spelling variations
            handle = Entrez.esearch(db="gene", term=gene_name + "[All Fields]")
            record = Entrez.read(handle)
            if record['IdList']:
                return True
            else:
                # If exact match and fuzzy match fail, try checking if it's a valid DNA sequence
                seq = Seq.Seq(gene_name)
                if seq.alphabet == Seq.Alphabet.IUPAC.unambiguous_dna:
                    return True
                else:
                    return False
    except Exception as e:
        print("Error occurred:", e)
        return False

xls = pd.ExcelFile("data/test.xlsx")
sheet_results = {}

# Iterate through each sheet and validate genes
for sheet_name in xls.sheet_names:
    # Read the data from the current sheet into a DataFrame
    df = pd.read_excel(xls, sheet_name=sheet_name)

    # Extract gene names from the DataFrame
    gene_names = df["Genes"].tolist()

    # Counter for the number of genes extracted
    total_genes = len(gene_names)

    # Counter for the number of validated genes
    validated_genes = 0

    # Validate each gene name
    for gene_name in gene_names:
        if is_valid_gene(gene_name):
            validated_genes += 1

    # Calculate the accuracy of valid genes
    if total_genes > 0:
        accuracy = (validated_genes / total_genes) * 100
    else:
        accuracy = 0

    # Store the results for the current sheet
    sheet_results[sheet_name] = {
        "total_genes": total_genes,
        "validated_genes": validated_genes,
        "accuracy": accuracy
    }

print("=============================================")

# Print the results for each sheet
for sheet_name, result in sheet_results.items():
    print(f"Sheet {sheet_name}:")
    print("Number of genes extracted:", result["total_genes"])
    print("Number of validated genes:", result["validated_genes"])
    print("Accuracy of valid genes:", result["accuracy"], "%")
    print()

print("=============================================")

# Initialize lists to store the data
sheets = []
total_genes_extracted = []
validated_genes_count = []
accuracy_percent = []

# Iterate through the sheet results dictionary
for sheet_name, result in sheet_results.items():
    # Append data to lists
    sheets.append(sheet_name)
    total_genes_extracted.append(result["total_genes"])
    validated_genes_count.append(result["validated_genes"])
    accuracy_percent.append(result["accuracy"])

# Create a dictionary with the data
data = {
    "Sheet": sheets,
    "Number of genes extracted": total_genes_extracted,
    "Number of validated genes": validated_genes_count,
    "Accuracy of valid genes (%)": accuracy_percent
}

# Create a DataFrame from the dictionary
df = pd.DataFrame(data)

# Check if path is correct and Excel file exists, adjust accordingly
file_path = "data/test_report.xlsx"

# Save the DataFrame to an existing Excel file, adding it as a new sheet
with pd.ExcelWriter(file_path, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
    df.to_excel(writer, sheet_name='Gene Validation Summary', index=False)

print(f"Gene validation summary has been appended to {file_path}.")