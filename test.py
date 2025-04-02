import pandas as pd

df = pd.read_csv("mito.csv")
print(df.head())  # Check the first few rows
print(df["Allele"].unique())  # Check unique allele formats
