from Bio import Entrez, Seq
import pandas as pd 

def is_valid_gene(gene_name):
    try:
        handle = Entrez.esearch(db="gene", term=gene_name)
        record = Entrez.read(handle)
        if record['IdList']:
            return True
        else:
            handle = Entrez.esearch(db="gene", term=gene_name + "[All Fields]")
            record = Entrez.read(handle)
            if record['IdList']:
                return True
            else:
                seq = Seq.Seq(gene_name)
                if seq.alphabet == Seq.Alphabet.IUPAC.unambiguous_dna:
                    return True
                else:
                    return False
    except Exception as e:
        print("Error occurred:", e)
        return False

def run_validation():
    # Set up Entrez
    Entrez.email = "your-email"

    xls = pd.ExcelFile("data/test.xlsx")
    sheet_results = {}

    for sheet_name in xls.sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet_name)
        gene_names = df["Genes"].tolist()
        total_genes = len(gene_names)
        validated_genes = 0
        for gene_name in gene_names:
            if is_valid_gene(gene_name):
                validated_genes += 1
        if total_genes > 0:
            accuracy = (validated_genes / total_genes) * 100
        else:
            accuracy = 0
        sheet_results[sheet_name] = {
            "total_genes": total_genes,
            "validated_genes": validated_genes,
            "accuracy": accuracy
        }

    create_summary(sheet_results)

def create_summary(sheet_results):
    sheets = []
    total_genes_extracted = []
    validated_genes_count = []
    accuracy_percent = []
    for sheet_name, result in sheet_results.items():
        sheets.append(sheet_name)
        total_genes_extracted.append(result["total_genes"])
        validated_genes_count.append(result["validated_genes"])
        accuracy_percent.append(result["accuracy"])
    data = {
        "Sheet": sheets,
        "Number of genes extracted": total_genes_extracted,
        "Number of validated genes": validated_genes_count,
        "Accuracy of valid genes (%)": accuracy_percent
    }
    df = pd.DataFrame(data)
    file_path = "data/test_report.xlsx"
    with pd.ExcelWriter(file_path, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
        df.to_excel(writer, sheet_name='Gene Validation Summary', index=False)
    print(f"Gene validation summary has been appended to {file_path}.")

if __name__ == "__main__":
    run_validation()