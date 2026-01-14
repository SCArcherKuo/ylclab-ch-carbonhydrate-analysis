# Module Organization

This document describes the structure of the carbohydrate analysis package.

## Project Structure

```
src/carbonhydrate_analysis/
├── __init__.py           # Package exports and version info
├── config.py             # Configuration constants
├── main.py               # Main entry point with examples
├── utils.py              # Utility functions
├── chebi_api.py          # ChEBI API client
├── pubchem_api.py        # PubChem API client
└── classification.py     # Carbohydrate classification logic
```

## Module Descriptions

### `main.py` - Main Entry Point
- Contains demonstration code showing how to use the package
- Provides example single and batch queries
- Helper functions for displaying results

**Run it:**
```bash
uv run python -m carbonhydrate_analysis.main
```

### `utils.py` - Utility Functions
**Purpose:** Helper functions for data processing and formatting

**Key Functions:**
- `extract_term_string()` - Extract strings from PubChem's StringWithMarkup format
- `extract_ontology_terms_from_node()` - Recursively extract ontology terms from classification nodes

### `chebi_api.py` - ChEBI API Client
**Purpose:** Interact with the ChEBI database API

**Key Classes:**
- `ChEBIClient` - Object-oriented client for ChEBI API with caching

**Key Functions:**
- `get_chebi_children()` - Get direct children of a ChEBI ID
- `get_main_groups()` - Get main groups (children with >1 children)

**Features:**
- Built-in caching to avoid repeated API calls
- Error handling and retry logic
- Clean separation of concerns

### `classification.py` - Carbohydrate Classification
**Purpose:** Classify compounds as carbohydrates based on ChEBI ontology

**Key Classes:**
- `CarbohydrateClassifier` - Implements 5-category classification system

**Key Functions:**
- `classify_carbohydrate()` - Main classification function

**Classification Categories:**
1. `main carbohydrate group` - Major carbohydrate groups (e.g., monosaccharide, oligosaccharide)
2. `other carbohydrate` - Other direct carbohydrate children
3. `main carbohydrate derivative group` - Major derivative groups (e.g., glycosides, amino sugars)
4. `other carbohydrate derivative` - Other derivative children
5. `other` - Under root but not carbohydrate or derivative
6. `None` - Not a carbohydrate

### `pubchem_api.py` - PubChem API Client
**Purpose:** Interact with the PubChem database API

**Key Classes:**
- `PubChemClient` - Object-oriented client for PubChem API

**Key Functions:**
- `get_compound_info_pubchem()` - Main API function (supports single and batch queries)

**Features:**
- Automatic identifier type detection (InChIKey vs SMILES)
- Batch processing with POST API for efficiency
- Rate limiting and error handling
- Two-step process:
  1. Resolve identifiers to CIDs
  2. Batch fetch properties, classifications, and synonyms

## Usage Examples

### Basic Usage

```python
from carbonhydrate_analysis import get_compound_info_pubchem

# Single query
result = get_compound_info_pubchem("RBNPOMFGQQGHHO-UHFFFAOYSA-N")
print(f"Name: {result['name']}")
print(f"Is Carbohydrate: {result['is_carbohydrate']}")

# Batch query
identifiers = [
    "RBNPOMFGQQGHHO-UHFFFAOYSA-N",
    "WQZGKKKJIJFFOK-GASJEMHNSA-N"
]
results = get_compound_info_pubchem(identifiers)
for result in results:
    if result:
        print(f"{result['name']}: {result['is_carbohydrate']}")
```

### Using Classes Directly

```python
from carbonhydrate_analysis import PubChemClient, ChEBIClient, CarbohydrateClassifier

# Create custom clients with different settings
pubchem = PubChemClient(timeout=60)
chebi = ChEBIClient()
classifier = CarbohydrateClassifier()

# Use clients
cid = pubchem.resolve_identifier_to_cid("RBNPOMFGQQGHHO-UHFFFAOYSA-N", "inchikey")
children = chebi.get_children(16646)  # carbohydrate ChEBI ID
main_class, subclass = classifier.classify(ontology_terms)
```

## Benefits of Reorganization

### 1. **Separation into Multiple Modules**
- **chebi_api.py** - ChEBI-specific functionality
- **pubchem_api.py** - PubChem-specific functionality
- **classification.py** - Classification logic
- **utils.py** - Shared utilities
- **main.py** - Entry point and examples

### 2. **Clear Code Organization**
- Each module has a clear purpose
- Related functions grouped together
- Section comments within modules
- Comprehensive docstrings

### 3. **Helper Functions in Utils**
- `extract_term_string()` - Data format conversion
- `extract_ontology_terms_from_node()` - Recursive tree traversal
- Reusable across modules

### 4. **Class-Based Structure**
- `ChEBIClient` - Encapsulates ChEBI API logic and caching
- `PubChemClient` - Encapsulates PubChem API logic
- `CarbohydrateClassifier` - Encapsulates classification rules
- Benefits:
  - State management (caching)
  - Configurable instances
  - Easier testing and mocking
  - Backward compatible with function-based API

## API Design

The package offers both **functional** and **object-oriented** interfaces:

**Functional (Simple):**
```python
from carbonhydrate_analysis import get_compound_info_pubchem
result = get_compound_info_pubchem("RBNPOMFGQQGHHO-UHFFFAOYSA-N")
```

**Object-Oriented (Advanced):**
```python
from carbonhydrate_analysis import PubChemClient
client = PubChemClient(timeout=60, rate_limit_delay=0.2)
result = client.get_compound_info_batch([752, 5793])
```

Both interfaces work seamlessly together!

## Testing

Run the main module to test functionality:

```bash
# Test both single and batch queries
uv run python -m carbonhydrate_analysis.main
```

## Future Improvements

Potential enhancements:
- Add async support for parallel API requests
- Implement result caching/database
- Add command-line interface (CLI)
- Create data processing pipeline for metabolomics datasets
- Add export functions (CSV, JSON, Excel)
