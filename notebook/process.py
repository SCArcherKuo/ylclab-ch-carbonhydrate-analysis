# %%
"""
Process metabolomics data and classify carbohydrates

This script provides functions to parse metabolomics Excel files and classify
compounds as carbohydrates using ChEBI ontology.
"""

# %%
import pandas as pd
from typing import List, Optional
from tqdm import tqdm
from carbonhydrate_analysis import get_compound_info_pubchem
from carbonhydrate_analysis.utils import extract_ontology_terms_from_node

# %%
def parse_metabolomics_xlsx(
    file_path: str,
    header_row: int = 4
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parse metabolomics Excel file and separate metadata and sample data columns.
    
    Parameters:
    -----------
    file_path : str
        Path to the Excel file containing metabolomics data
    header_row : int, optional
        Row index where the header is located (default: 4)
    
    Returns:
    --------
    tuple of (full_df, metadata_df, sample_data_df)
        full_df: Complete dataframe with all columns
        metadata_df: Dataframe with only metadata columns
        sample_data_df: Dataframe with sample data columns (those containing 'CH_Algae')
    
    Examples:
    ---------
    >>> full_df, metadata_df, sample_df = parse_metabolomics_xlsx(
    ...     'data/raw/20251126_CH_Algae_T3_metabolomes_NEG_Area.xlsx'
    ... )
    >>> print(f"Total columns: {len(full_df.columns)}")
    >>> print(f"Metadata columns: {len(metadata_df.columns)}")
    >>> print(f"Sample columns: {len(sample_df.columns)}")
    """
    # Read the Excel file with proper header
    df = pd.read_excel(file_path, header=header_row)
    
    # Identify sample columns (those containing 'CH_Algae' in the name)
    sample_columns = [col for col in df.columns if 'CH_Algae' in str(col)]
    
    # Identify metadata columns (all columns except sample columns)
    metadata_columns = [col for col in df.columns if col not in sample_columns]
    
    # Create separate dataframes
    metadata_df = df[metadata_columns].copy()
    sample_data_df = df[sample_columns].copy()
    
    return df, metadata_df, sample_data_df

# %%
def classify_carbohydrates(
    metadata_df: pd.DataFrame,
    inchikey_column: str = 'INCHIKEY'
) -> pd.DataFrame:
    """
    Classify compounds as carbohydrates using InChIKey and PubChem API.
    
    This function takes a metadata dataframe and adds three new columns:
    - 'Is Carbohydrate': Boolean indicating if the compound is a carbohydrate
    - 'Main Class': Main classification category
    - 'Subclass': Specific subclass name
    
    Parameters:
    -----------
    metadata_df : pd.DataFrame
        Dataframe containing metabolite metadata with an INCHIKEY column
    inchikey_column : str, optional
        Name of the column containing InChIKey information (default: 'INCHIKEY')
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with three additional columns: 'Is Carbohydrate', 'Main Class', 'Subclass'
    
    Examples:
    ---------
    >>> metadata_df = parse_metabolomics_xlsx('data.xlsx')[1]
    >>> classified_df = classify_carbohydrates(metadata_df)
    >>> print(f"Carbohydrates found: {classified_df['Is Carbohydrate'].sum()}")
    """
    # Create a copy to avoid modifying the original
    result_df = metadata_df.copy()
    
    # Initialize new columns
    result_df['Is Carbohydrate'] = False
    result_df['Main Class'] = None
    result_df['Subclass'] = None
    
    # Collect all valid InChIKeys
    valid_inchikeys = []
    inchikey_to_idx = {}
    
    for idx, row in result_df.iterrows():
        inchikey_value = row.get(inchikey_column)
        
        # Skip if InChIKey is missing or NaN
        if pd.notna(inchikey_value) and isinstance(inchikey_value, str) and inchikey_value.strip():
            inchikey = inchikey_value.strip()
            valid_inchikeys.append(inchikey)
            if inchikey not in inchikey_to_idx:
                inchikey_to_idx[inchikey] = []
            inchikey_to_idx[inchikey].append(idx)
    
    if not valid_inchikeys:
        print("No valid InChIKeys found in the dataset")
        return result_df
    
    print(f"Found {len(valid_inchikeys)} valid InChIKeys ({len(inchikey_to_idx)} unique)")
    
    # Get unique InChIKeys for batch query
    unique_inchikeys = list(inchikey_to_idx.keys())
    
    # Query PubChem in batch
    print(f"\nQuerying PubChem for {len(unique_inchikeys)} unique compounds...")
    pubchem_results = get_compound_info_pubchem(unique_inchikeys, identifier_type='inchikey')
    
    # Process results with progress bar
    print("\nProcessing results...")
    for inchikey, compound_info in tqdm(
        zip(unique_inchikeys, pubchem_results),
        total=len(unique_inchikeys),
        desc="Classifying compounds"
    ):
        # Get all rows with this InChIKey
        row_indices = inchikey_to_idx[inchikey]
        
        if compound_info is not None:
            is_carb = compound_info.get('is_carbohydrate', False)
            main_class = compound_info.get('carbohydrate_main_class')
            subclass = compound_info.get('carbohydrate_subclass')
            
            # Update all rows with this InChIKey
            for idx in row_indices:
                result_df.loc[idx, 'Is Carbohydrate'] = is_carb  # type: ignore
                if is_carb:
                    result_df.loc[idx, 'Main Class'] = main_class  # type: ignore
                    result_df.loc[idx, 'Subclass'] = subclass  # type: ignore
    
    return result_df

# %%
def process_metabolomics_file(
    file_path: str,
    header_row: int = 4,
    inchikey_column: str = 'INCHIKEY'
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Complete processing pipeline: parse file and classify carbohydrates.
    
    This convenience function combines parsing and classification into one step.
    
    Parameters:
    -----------
    file_path : str
        Path to the Excel file containing metabolomics data
    header_row : int, optional
        Row index where the header is located (default: 4)
    inchikey_column : str, optional
        Name of the column containing InChIKey information (default: 'INCHIKEY')
    
    Returns:
    --------
    tuple of (full_df, classified_metadata_df, sample_data_df)
        full_df: Complete dataframe with all columns
        classified_metadata_df: Metadata with carbohydrate classification columns
        sample_data_df: Sample data columns
    
    Examples:
    ---------
    >>> full_df, metadata_df, sample_df = process_metabolomics_file(
    ...     'data/raw/20251126_CH_Algae_T3_metabolomes_NEG_Area.xlsx'
    ... )
    >>> carb_count = metadata_df['Is Carbohydrate'].sum()
    >>> print(f"Found {carb_count} carbohydrates")
    """
    # Parse the file
    full_df, metadata_df, sample_data_df = parse_metabolomics_xlsx(file_path, header_row)
    
    # Classify carbohydrates using InChIKey
    classified_metadata_df = classify_carbohydrates(metadata_df, inchikey_column)
    
    return full_df, classified_metadata_df, sample_data_df

# %%
# Example usage
if __name__ == '__main__':
    import os
    import glob
    
    # Determine the correct path based on where the script is run from
    if os.path.exists('../data/raw'):
        raw_dir = '../data/raw'
        output_dir = '../data/processed'
    else:
        raw_dir = 'data/raw'
        output_dir = 'data/processed'
    
    # Find all xlsx files in data/raw
    xlsx_files = glob.glob(os.path.join(raw_dir, '*.xlsx'))
    
    if not xlsx_files:
        print(f"No Excel files found in {raw_dir}")
        exit(1)
    
    print(f"Found {len(xlsx_files)} Excel file(s) to process:\n")
    for i, file in enumerate(xlsx_files, 1):
        print(f"  {i}. {os.path.basename(file)}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each file
    total_carbohydrates = 0
    for file_idx, file_path in enumerate(xlsx_files, 1):
        print(f"\n{'='*80}")
        print(f"Processing file {file_idx}/{len(xlsx_files)}: {os.path.basename(file_path)}")
        print('='*80)
        
        try:
            full_df, metadata_df, sample_df = process_metabolomics_file(file_path)
            
            print(f"\nFile structure:")
            print(f"  Total columns: {len(full_df.columns)}")
            print(f"  Metadata columns: {len(metadata_df.columns)}")
            print(f"  Sample columns: {len(sample_df.columns)}")
            
            carb_count = metadata_df['Is Carbohydrate'].sum()
            print(f"\nMetabolite classification:")
            print(f"  Total metabolites: {len(metadata_df)}")
            print(f"  Carbohydrates found: {carb_count}")
            total_carbohydrates += carb_count
            
            if carb_count > 0:
                print(f"\nCarbohydrate breakdown:")
                class_counts = metadata_df[metadata_df['Is Carbohydrate']]['Main Class'].value_counts()
                for class_name, count in class_counts.items():
                    print(f"  {class_name}: {count}")
            
            # Extract base filename without extension
            base_filename = os.path.splitext(os.path.basename(file_path))[0]
            
            # Save classified metadata
            metadata_output = os.path.join(output_dir, f"{base_filename}_classified_metadata.csv")
            metadata_df.to_csv(metadata_output, index=False)
            print(f"\n✓ Saved classified metadata to: {metadata_output}")
            
            # Save sample data
            sample_output = os.path.join(output_dir, f"{base_filename}_sample_data.csv")
            sample_df.to_csv(sample_output, index=False)
            print(f"✓ Saved sample data to: {sample_output}")
            
            # Save carbohydrates only (if any found)
            if carb_count > 0:
                carb_only = metadata_df[metadata_df['Is Carbohydrate']].copy()
                carb_output = os.path.join(output_dir, f"{base_filename}_carbohydrates_only.csv")
                carb_only.to_csv(carb_output, index=False)
                print(f"✓ Saved carbohydrates only to: {carb_output}")
            
        except Exception as e:
            print(f"✗ Error processing {os.path.basename(file_path)}: {str(e)}")
            continue
    
    print(f"\n{'='*80}")
    print(f"SUMMARY")
    print('='*80)
    print(f"Files processed: {len(xlsx_files)}")
    print(f"Total carbohydrates found: {total_carbohydrates}")
    print(f"All results saved to {output_dir}/")

# %%
