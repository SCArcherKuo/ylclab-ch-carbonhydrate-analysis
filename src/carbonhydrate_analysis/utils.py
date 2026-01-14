"""
Utility Functions

This module provides helper functions for data processing and formatting.
"""

from typing import Any


def extract_term_string(term: Any) -> str:
    """
    Extract string from StringWithMarkup format or return as-is.
    
    This function handles PubChem's StringWithMarkup format which is used
    in classification hierarchies and ontology terms.
    
    Parameters:
    -----------
    term : any
        Term that might be in StringWithMarkup format or plain string
    
    Returns:
    --------
    str
        Extracted string value
    
    Examples:
    ---------
    >>> extract_term_string("simple string")
    'simple string'
    >>> extract_term_string({'StringWithMarkup': {'String': 'formatted'}})
    'formatted'
    """
    if isinstance(term, dict):
        # Check for StringWithMarkup format
        if 'StringWithMarkup' in term:
            string_data = term['StringWithMarkup']
            if isinstance(string_data, dict) and 'String' in string_data:
                return string_data['String']
        # Try direct String key
        if 'String' in term:
            return term['String']
    return str(term)


def extract_ontology_terms_from_node(node: dict, terms_list: list) -> None:
    """
    Recursively extract ontology terms from classification nodes.
    
    This function traverses the hierarchical structure of PubChem
    classification data and extracts all ontology term names.
    
    Parameters:
    -----------
    node : dict
        Classification node from PubChem hierarchy
    terms_list : list
        List to append extracted terms to (modified in-place)
    
    Examples:
    ---------
    >>> terms = []
    >>> node = {'Information': {'Name': 'carbohydrate'}}
    >>> extract_ontology_terms_from_node(node, terms)
    >>> print(terms)
    ['carbohydrate']
    """
    if not isinstance(node, dict):
        return
    
    # Extract name from current node
    if 'Information' in node:
        info = node['Information']
        if 'Name' in info:
            terms_list.append(info['Name'])
    
    # Recursively process children
    if 'Children' in node and 'Node' in node['Children']:
        children = node['Children']['Node']
        if isinstance(children, list):
            for child in children:
                extract_ontology_terms_from_node(child, terms_list)
        else:
            extract_ontology_terms_from_node(children, terms_list)
