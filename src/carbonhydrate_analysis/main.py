"""
Carbohydrate Filtering for Metabolomics Data

Main entry point and demonstration of the carbohydrate analysis functionality.
"""

from typing import Dict, Any
from .pubchem_api import get_compound_info_pubchem


def main():
    """
    Main function for testing the compound lookup functionality.
    
    Demonstrates both single and batch query capabilities with example compounds.
    """
    print("Carbohydrate Filtering Tool")
    print("=" * 50)
    
    # ==========================================================================
    # Test 1: Single InChIKey Query
    # ==========================================================================
    
    test_inchikey = "RBNPOMFGQQGHHO-UHFFFAOYSA-N"
    print(f"\n[Test 1] Single InChIKey query: {test_inchikey}")
    
    result = get_compound_info_pubchem(test_inchikey)
    
    if result and isinstance(result, dict):
        _display_compound_result(result)
    else:
        print("\n✗ Failed to retrieve compound information")
    
    # ==========================================================================
    # Test 2: Batch InChIKey Query
    # ==========================================================================
    
    print("\n" + "=" * 50)
    test_inchikeys = [
        "RBNPOMFGQQGHHO-UHFFFAOYSA-N",  # Glyceric acid
        "WQZGKKKJIJFFOK-GASJEMHNSA-N"   # Glucose
    ]
    print(f"\n[Test 2] Batch InChIKey query: {len(test_inchikeys)} identifiers")
    
    results = get_compound_info_pubchem(test_inchikeys, identifier_type='inchikey')
    
    if isinstance(results, list):
        _display_batch_results(results, test_inchikeys)


def _display_compound_result(result: Dict[str, Any]) -> None:
    """
    Display detailed information for a single compound result.
    
    Parameters:
    -----------
    result : dict
        Compound information dictionary
    """
    print(f"\n✓ Compound found:")
    print(f"  Name: {result['name']}")
    print(f"  Formula: {result['formula']}")
    print(f"  PubChem CID: {result['pubchem_cid']}")
    print(f"  Is Carbohydrate: {result['is_carbohydrate']}")
    
    if result['is_carbohydrate']:
        print(f"  Main Class: {result['carbohydrate_main_class']}")
        print(f"  Subclass: {result['carbohydrate_subclass']}")
    
    if result['chebi_ontology']:
        print(f"\n  ChEBI Ontology Terms ({len(result['chebi_ontology'])}):")
        for term in result['chebi_ontology'][:10]:
            print(f"    - {term}")
        if len(result['chebi_ontology']) > 10:
            print(f"    ... and {len(result['chebi_ontology']) - 10} more")


def _display_batch_results(results: list, identifiers: list) -> None:
    """
    Display summary information for batch query results.
    
    Parameters:
    -----------
    results : list
        List of compound information dictionaries
    identifiers : list
        Original list of identifiers queried
    """
    successful = len([r for r in results if r])
    print(f"\n✓ Retrieved {successful} out of {len(identifiers)} compounds")
    
    for i, res in enumerate(results):
        if res:
            print(f"\n  [{i+1}] {res['name']}")
            print(f"      CID: {res['pubchem_cid']}, Is Carbohydrate: {res['is_carbohydrate']}")
            if res['is_carbohydrate']:
                print(f"      Class: {res['carbohydrate_main_class']} - {res['carbohydrate_subclass']}")
        else:
            print(f"\n  [{i+1}] Failed to retrieve")


if __name__ == "__main__":
    main()
