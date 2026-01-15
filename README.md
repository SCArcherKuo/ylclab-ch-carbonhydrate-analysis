# Carbohydrate Analysis for Metabolomics Data

A high-performance Python package for identifying and classifying carbohydrates in metabolomics data using PubChem and ChEBI ontology databases.

[![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Features

- **High-Performance Classification** - Optimized workflow with 99% bandwidth reduction and 50% faster processing
- **Dual Classification Methods** - Efficient ChEBI ancestry-based method with automatic fallback to PubChem
- **Intelligent Caching** - Persistent cache system for ChEBI ontology relationships
- **Batch Processing** - Efficiently process multiple compounds simultaneously
- **Five-Category Classification System** - Detailed carbohydrate categorization
- **Automatic InChIKey Resolution** - Seamlessly converts InChIKeys to compound information

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/SCArcherKuo/ylclab-ch-carbonhydrate-analysis.git
cd ylclab-ch-carbonhydrate-analysis

# Install dependencies using uv (recommended)
uv sync

# Or using pip
pip install -e .
```

### Basic Usage

```python
from carbonhydrate_analysis import get_compound_info_pubchem

# Single compound query
result = get_compound_info_pubchem("WQZGKKKJIJFFOK-GASJEMHNSA-N")

print(f"Compound: {result['name']}")
print(f"Is Carbohydrate: {result['is_carbohydrate']}")
print(f"Main Class: {result['carbohydrate_main_class']}")
print(f"Subclass: {result['carbohydrate_subclass']}")

# Batch processing
inchikeys = [
    "WQZGKKKJIJFFOK-GASJEMHNSA-N",  # D-Glucose
    "RBNPOMFGQQGHHO-UHFFFAOYSA-N"   # Glyceric acid
]
results = get_compound_info_pubchem(inchikeys, identifier_type='inchikey')

for result in results:
    if result and result['is_carbohydrate']:
        print(f"{result['name']}: {result['carbohydrate_subclass']}")
```

### Process Metabolomics Files

```python
from notebook.process import process_metabolomics_file

# Process an MS-DIAL output stored in a Excel file
full_df, metadata_df, sample_df = process_metabolomics_file(
    'data/raw/your_metabolomics_data.xlsx',
    header_row=4,
    inchikey_column='INCHIKEY'
)

# Count carbohydrates
carb_count = metadata_df['Is Carbohydrate'].sum()
print(f"Found {carb_count} carbohydrates")

# Filter by class
monosaccharides = metadata_df[metadata_df['Subclass'] == 'monosaccharide']
```

## Classification System

The package uses a five-category classification system based on ChEBI ontology (CHEBI:78616):

### Categories

1. **Main Carbohydrate Group** - Major carbohydrate classes with significant ontology branching
   - Examples: monosaccharide, disaccharide, polysaccharide, oligosaccharide

2. **Other Carbohydrate** - Carbohydrates not in main groups
   - Direct children of CHEBI:16646 with fewer branches

3. **Main Carbohydrate Derivative Group** - Major derivative classes
   - Examples: amino sugar, deoxy sugar, sugar phosphate

4. **Other Carbohydrate Derivative** - Derivatives not in main groups
   - Direct children of CHEBI:63299 with fewer branches

5. **Other** - Under carbohydrates root but not classified as carbohydrate or derivative

### Classification Logic

```
carbohydrates and carbohydrate derivatives (CHEBI:78616)
├── carbohydrate (CHEBI:16646)
│   ├── Main Groups (>1 children)
│   │   ├── monosaccharide
│   │   ├── disaccharide
│   │   └── polysaccharide
│   └── Other Groups (≤1 children)
└── carbohydrate derivative (CHEBI:63299)
    ├── Main Groups (>1 children)
    │   ├── amino sugar
    │   ├── deoxy sugar
    │   └── sugar phosphate
    └── Other Groups (≤1 children)
```

## Architecture

### Optimized Workflow

The package uses an intelligent dual-method approach:

**Fast Path (ChEBI Ancestry - 80%+ of compounds):**
```
InChIKey → CID → Synonyms → Extract ChEBI ID → Query ChEBI Ancestry → Classify
```
- 99% less bandwidth (KB vs GB)
- 50% faster
- Highly cacheable

**Fallback Path (PubChem Classification):**
```
InChIKey → CID → Full Classification → Extract ChEBI Terms → Classify
```
- Comprehensive coverage
- No functionality lost
- Same accuracy

### Module Structure

```
src/carbonhydrate_analysis/
├── __init__.py              # Package initialization
├── pubchem_api.py          # PubChem API client with optimization
├── chebi_api.py            # ChEBI API client with ancestry support
├── classification.py       # Classification logic (both methods)
├── utils.py                # Utility functions
├── config.py               # Configuration constants
├── main.py                 # Demo script
└── cache_stats.py          # Cache statistics utility
```

## API Reference

### Main Functions

#### `get_compound_info_pubchem(identifier, identifier_type='auto')`

Fetch compound information from PubChem using InChIKey or SMILES.

**Parameters:**
- `identifier` (str or list): Single identifier or list of identifiers
- `identifier_type` (str): 'inchikey', 'smiles', or 'auto' (default)

**Returns:**
- Dict or list of dicts with compound information:
  - `pubchem_cid`: PubChem Compound ID
  - `name`: IUPAC name
  - `formula`: Molecular formula
  - `molecular_weight`: Molecular weight
  - `inchi`: InChI string
  - `inchikey`: InChIKey
  - `smiles`: Canonical SMILES
  - `chebi_id`: ChEBI ID (if found)
  - `is_carbohydrate`: Boolean flag
  - `carbohydrate_main_class`: Main classification category
  - `carbohydrate_subclass`: Specific subclass name
  - `classifications`: Classification hierarchies (empty if ChEBI method used)
  - `chebi_ontology`: ChEBI ontology terms (empty if ChEBI method used)

### ChEBI API Functions

#### `get_all_ancestors(chebi_id)`

Get all ancestor ChEBI IDs by recursively traversing parent relationships.

**Parameters:**
- `chebi_id` (int): ChEBI ID (without 'CHEBI:' prefix)

**Returns:**
- List of int: All ancestor ChEBI IDs (includes the starting ID)

#### `get_chebi_children(chebi_id)`

Get direct children of a ChEBI ID.

**Parameters:**
- `chebi_id` (int): ChEBI ID (without 'CHEBI:' prefix)

**Returns:**
- List of dict: Direct children with 'is a' relationship

### Classification Functions

#### `classify_by_chebi_ancestry(chebi_id)`

Classify a compound based on its ChEBI ID ancestry (efficient method).

**Parameters:**
- `chebi_id` (int): ChEBI ID (without 'CHEBI:' prefix)

**Returns:**
- Tuple of (main_class, subclass)

#### `classify_carbohydrate(chebi_ontology)`

Classify a compound based on its ChEBI ontology terms (fallback method).

**Parameters:**
- `chebi_ontology` (list): List of ChEBI ontology terms

**Returns:**
- Tuple of (main_class, subclass)


## Testing

### Run Demo Script

```bash
# Test with example compounds
python -m carbonhydrate_analysis.main
```

### Process Your Data

```bash
# Process Excel files in data/raw/
uv run notebook/process.py
```

Expected output:
```
Found 2 Excel file(s) to process:
  1. your_data.xlsx

Processing file 1/1: your_data.xlsx
==========================================
Resolving 150 inchikeys to CIDs...
Found 145 valid CIDs, fetching properties...
Processing compounds: 100%|██████████| 145/145

Metabolite classification:
  Total metabolites: 150
  Carbohydrates found: 42

Carbohydrate breakdown:
  main carbohydrate group: 35
  other carbohydrate: 5
  main carbohydrate derivative group: 2
```


## Requirements

- Python 3.13+
- pandas
- requests
- tqdm
- openpyxl (for Excel file processing)
- numpy
- pyasn1 (for future ASNT format support)


## Data Sources

- **PubChem**: Chemical structure database (https://pubchem.ncbi.nlm.nih.gov/)
- **ChEBI**: Chemical Entities of Biological Interest (https://www.ebi.ac.uk/chebi/)
