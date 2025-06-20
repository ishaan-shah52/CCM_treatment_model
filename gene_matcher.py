import pandas as pd

def extract_and_match_columns(file1, sheet1, column1, output_csv):
    # Load the Excel files and extract the specific columns
    my_network = pd.read_excel(file1, sheet_name=sheet1, skiprows=[0])
    all_drugs_CSV = pd.read_csv('alldrugbank.csv') # from DrugBank website
    
    # Extract the columns
    col1_data = my_network[column1].astype(str)  # Convert to string for matching
    col2_data = all_drugs_CSV['Gene Name'].astype(str)  # Convert to string for matching
    
    # Create a list to store the results
    results = []
    
    # Iterate through each value in the first column
    for value in col1_data:
        # Check if the value is present anywhere in the second column
        match = col2_data[col2_data.str.contains(value, case=False, na=False)]
        if not match.empty:
            results.append(value)
        else:
            results.append("None")
    
    # Create a DataFrame for the results
    results_df = pd.DataFrame({column1: col1_data, 'Match': results})
    
    # Save the results to a CSV file
    results_df.to_csv(output_csv, index=False)

# Example usage
extract_and_match_columns('CCMNetIS.xlsx', 'species', 'ID', 'matched.csv')
