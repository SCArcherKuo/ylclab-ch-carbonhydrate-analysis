"""
Carbohydrate Classification

This module provides functions to classify compounds as carbohydrates based on
their ChEBI ontology terms and hierarchical relationships.
"""

from typing import List, Any, Tuple, Optional
from .chebi_api import get_chebi_children, get_main_groups
from .utils import extract_term_string
from . import config


class CarbohydrateClassifier:
    """
    Classifier for identifying and categorizing carbohydrates.
    
    This class implements the classification logic for determining whether
    a compound is a carbohydrate and categorizing it into one of five main
    categories based on ChEBI ontology.
    """
    
    def __init__(self):
        """Initialize the carbohydrate classifier."""
        self.root_id = config.CHEBI_ROOT_ID
        self.carbohydrate_id = config.CHEBI_CARBOHYDRATE_ID
        self.carbohydrate_derivative_id = config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID
    
    def classify(self, chebi_ontology: List[Any]) -> Tuple[Optional[str], Optional[str]]:
        """
        Classify a compound based on its ChEBI ontology terms.
        
        Parameters:
        -----------
        chebi_ontology : list
            List of ChEBI ontology terms (may be strings or StringWithMarkup dicts)
        
        Returns:
        --------
        tuple of (main_class, subclass)
            main_class: One of ['main carbohydrate group', 'other carbohydrate', 
                               'main carbohydrate derivative group', 'other carbohydrate derivative', 
                               'other', None]
            subclass: Specific subclass name or None
        
        Examples:
        ---------
        >>> classifier = CarbohydrateClassifier()
        >>> ontology = ['carbohydrates and carbohydrate derivatives', 'carbohydrate', 'monosaccharide']
        >>> main_class, subclass = classifier.classify(ontology)
        >>> print(main_class, subclass)
        'main carbohydrate group' 'monosaccharide'
        """
        if not chebi_ontology:
            return None, None
        
        # Extract clean strings from all terms
        clean_terms = [extract_term_string(term) for term in chebi_ontology]
        clean_terms_lower = [t.lower() for t in clean_terms]
        
        # Check if compound belongs to carbohydrates and carbohydrate derivatives (CHEBI:78616)
        if 'carbohydrates and carbohydrate derivatives' not in clean_terms_lower:
            return None, None
        
        # Get main groups for carbohydrate and carbohydrate derivative
        carb_main_groups = get_main_groups(self.carbohydrate_id)
        carb_deriv_main_groups = get_main_groups(self.carbohydrate_derivative_id)
        
        carb_main_groups_lower = [g.lower() for g in carb_main_groups]
        carb_deriv_main_groups_lower = [g.lower() for g in carb_deriv_main_groups]
        
        # Check for carbohydrate (CHEBI:16646)
        if 'carbohydrate' in clean_terms_lower:
            return self._classify_carbohydrate(
                clean_terms, clean_terms_lower,
                carb_main_groups, carb_main_groups_lower
            )
        
        # Check for carbohydrate derivative (CHEBI:63299)
        if 'carbohydrate derivative' in clean_terms_lower:
            return self._classify_carbohydrate_derivative(
                clean_terms, clean_terms_lower,
                carb_deriv_main_groups, carb_deriv_main_groups_lower
            )
        
        # If under root but not in carbohydrate or carbohydrate derivative
        return self._classify_other(clean_terms, clean_terms_lower)
    
    def _classify_carbohydrate(
        self,
        clean_terms: List[str],
        clean_terms_lower: List[str],
        main_groups: List[str],
        main_groups_lower: List[str]
    ) -> Tuple[str, str]:
        """Classify carbohydrate into main or other group."""
        # Check if matches main carbohydrate groups
        for i, term_lower in enumerate(clean_terms_lower):
            if term_lower in main_groups_lower:
                idx = main_groups_lower.index(term_lower)
                return 'main carbohydrate group', main_groups[idx]
        
        # If no main group match, it's "other carbohydrate"
        # Get all direct children and find which one matches
        carb_children = get_chebi_children(self.carbohydrate_id)
        carb_children_names = [child['init_name'].lower() for child in carb_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in carb_children_names:
                # Return the original capitalization from the ontology
                return 'other carbohydrate', term
        
        return 'other carbohydrate', 'carbohydrate'
    
    def _classify_carbohydrate_derivative(
        self,
        clean_terms: List[str],
        clean_terms_lower: List[str],
        main_groups: List[str],
        main_groups_lower: List[str]
    ) -> Tuple[str, str]:
        """Classify carbohydrate derivative into main or other group."""
        # Check if matches main carbohydrate derivative groups
        for i, term_lower in enumerate(clean_terms_lower):
            if term_lower in main_groups_lower:
                idx = main_groups_lower.index(term_lower)
                return 'main carbohydrate derivative group', main_groups[idx]
        
        # If no main group match, it's "other carbohydrate derivative"
        # Get all direct children and find which one matches
        carb_deriv_children = get_chebi_children(self.carbohydrate_derivative_id)
        carb_deriv_children_names = [child['init_name'].lower() for child in carb_deriv_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in carb_deriv_children_names:
                # Return the original capitalization from the ontology
                return 'other carbohydrate derivative', term
        
        return 'other carbohydrate derivative', 'carbohydrate derivative'
    
    def _classify_other(
        self,
        clean_terms: List[str],
        clean_terms_lower: List[str]
    ) -> Tuple[str, str]:
        """Classify as 'other' category."""
        # Get all direct children and find which one matches
        root_children = get_chebi_children(self.root_id)
        root_children_names = [child['init_name'].lower() for child in root_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in root_children_names and term_lower != 'carbohydrates and carbohydrate derivatives':
                # Return the original capitalization from the ontology
                return 'other', term
        
        return 'other', 'carbohydrates and carbohydrate derivatives'


# =============================================================================
# Module-level convenience function for backward compatibility
# =============================================================================

# Global classifier instance
_default_classifier = CarbohydrateClassifier()


def classify_carbohydrate(chebi_ontology: List[Any]) -> Tuple[Optional[str], Optional[str]]:
    """
    Classify a compound based on its ChEBI ontology terms.
    
    Convenience function that uses the default CarbohydrateClassifier.
    
    Parameters:
    -----------
    chebi_ontology : list
        List of ChEBI ontology terms (may be strings or StringWithMarkup dicts)
    
    Returns:
    --------
    tuple of (main_class, subclass)
        main_class: One of ['main carbohydrate group', 'other carbohydrate', 
                           'main carbohydrate derivative group', 'other carbohydrate derivative', 
                           'other', None]
        subclass: Specific subclass name or None
    """
    return _default_classifier.classify(chebi_ontology)
