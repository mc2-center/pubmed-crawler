import pandas as pd

# Load the first CSV file into a DataFrame (NEW MANIFEST)
file1 = './csvs/feb/2025-02-04_publications-manifest-nodec.csv'
df1 = pd.read_csv(file1)

# Load the second CSV file into another DataFrame
file2 = './csvs/2025-01-06_publications-manifest.csv'
df2 = pd.read_csv(file2)

# Identify "Pubmed Id" values from df2 that need to be removed from df1
pubmed_ids_to_remove = set(df2['Pubmed Id'])
# Remove rows from df1 that have matching "Pubmed Id" values in df2
df1 = df1[~df1['Pubmed Id'].isin(pubmed_ids_to_remove)]

# Save the modified DataFrame back to a new CSV file
output_file = './csvs/feb/2025-02-04_publications-manifest-final.csv'
df1.to_csv(output_file, index=False)

print(f"Rows with matching Pubmed Id removed. Saved to {output_file}")