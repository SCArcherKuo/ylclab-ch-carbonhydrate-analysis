"""
PubChem API Client

This module provides functions for interacting with the PubChem API to retrieve
compound information, properties, classifications, and synonyms.
"""

import requests
import json
import time
import re
from typing import List, Dict, Any, Union, Optional
from tqdm import tqdm
from . import config
from .classification import classify_carbohydrate, classify_by_chebi_ancestry
from .utils import extract_ontology_terms_from_node


class PubChemClient:
    """
    Client for interacting with PubChem API.
    
    This class provides methods to query PubChem database for compound
    information using various identifiers (InChIKey, SMILES, CID).
    """
    
    def __init__(self, base_url: Optional[str] = None, timeout: Optional[int] = None, rate_limit_delay: Optional[float] = None):
        """
        Initialize PubChem API client.
        
        Parameters:
        -----------
        base_url : str, optional
            Base URL for PubChem API (defaults to config.PUBCHEM_BASE_URL)
        timeout : int, optional
            Request timeout in seconds (defaults to config.API_TIMEOUT)
        rate_limit_delay : float, optional
            Delay between requests in seconds (defaults to config.RATE_LIMIT_DELAY)
        """
        self.base_url = base_url or config.PUBCHEM_BASE_URL
        self.timeout = timeout or config.API_TIMEOUT
        self.rate_limit_delay = rate_limit_delay or config.RATE_LIMIT_DELAY
        self.max_synonyms = config.MAX_SYNONYMS
        self.chunk_size = 512
    
    # =========================================================================
    # Identifier Resolution
    # =========================================================================
    
    def resolve_identifier_to_cid(self, identifier: str, identifier_type: str) -> Optional[int]:
        """
        Resolve an identifier (InChIKey or SMILES) to a PubChem CID.
        
        Parameters:
        -----------
        identifier : str
            The identifier string
        identifier_type : str
            Type of identifier: 'inchikey' or 'smiles'
        
        Returns:
        --------
        int or None
            PubChem CID if found, None otherwise
        """
        search_url = f"{self.base_url}/compound/{identifier_type}/cids/JSON"
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        payload = f"{identifier_type}={identifier}"
        
        try:
            response = requests.post(search_url, data=payload, headers=headers, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                    cids = data['IdentifierList']['CID']
                    if cids:
                        return cids[0]  # Take first CID
            return None
        except Exception as e:
            print(f"Error resolving identifier: {str(e)}")
            return None
    
    def resolve_identifiers_to_cids(
        self,
        identifiers: List[str],
        identifier_type: str
    ) -> Dict[str, Optional[int]]:
        """
        Resolve multiple identifiers to CIDs.
        
        Parameters:
        -----------
        identifiers : list of str
            List of identifiers
        identifier_type : str
            Type of identifier: 'inchikey' or 'smiles'
        
        Returns:
        --------
        dict
            Mapping of identifier to CID (or None if not found)
        """
        print(f"Resolving {len(identifiers)} {identifier_type}s to CIDs...")
        
        identifier_to_cid = {}
        
        for identifier in identifiers:
            cid = self.resolve_identifier_to_cid(identifier, identifier_type)
            identifier_to_cid[identifier] = cid
            
            if cid:
                print(f"  {identifier[:30]}... -> CID {cid}")
            else:
                print(f"  {identifier[:30]}... -> No CID found")
            
            # Small delay to avoid rate limiting
            time.sleep(0.1)
        
        return identifier_to_cid
    
    # =========================================================================
    # Batch Properties Fetching
    # =========================================================================
    
    def get_properties(self, cids: List[int]) -> Dict[int, Dict[str, Any]]:
        """
        Get properties for multiple CIDs using POST API.
        
        Parameters:
        -----------
        cids : list of int
            List of PubChem CIDs
        
        Returns:
        --------
        dict
            Mapping of CID to properties dictionary
        """
        properties = {}
        
        for i in range(0, len(cids), self.chunk_size):
            chunk_cids = cids[i:i + self.chunk_size]
            cid_list = ','.join(map(str, chunk_cids))
            
            try:
                props_url = f"{self.base_url}/compound/cid/property/MolecularFormula,MolecularWeight,InChI,InChIKey,SMILES,IUPACName/JSON"
                props_response = requests.post(
                    props_url,
                    data=f"cid={cid_list}",
                    headers={'Content-Type': 'application/x-www-form-urlencoded'},
                    timeout=self.timeout
                )
                
                if props_response.status_code == 200:
                    props_data = props_response.json()
                    if 'PropertyTable' in props_data and 'Properties' in props_data['PropertyTable']:
                        for props in props_data['PropertyTable']['Properties']:
                            cid = props.get('CID')
                            if cid:
                                properties[cid] = props
            except Exception as e:
                print(f"Error fetching properties: {str(e)}")
        
        return properties
    
    def get_classification(self, cid: int) -> List[Dict[str, Any]]:
        """
        Get classification hierarchies for a CID.
        
        Parameters:
        -----------
        cid : int
            PubChem CID
        
        Returns:
        --------
        list
            List of classification hierarchies
        """
        class_url = f"{self.base_url}/compound/cid/{cid}/classification/JSON"
        
        try:
            response = requests.get(class_url, timeout=self.timeout)
            
            if response.status_code == 200:
                class_data = response.json()
                if 'Hierarchies' in class_data and 'Hierarchy' in class_data['Hierarchies']:
                    return [h for h in class_data['Hierarchies']['Hierarchy'] if 'Node' in h]
        except json.JSONDecodeError:
            print(f"Warning: Classification data has JSON errors for CID {cid}, skipping...")
        except Exception as e:
            print(f"Error fetching classification: {str(e)}")
        
        return []
    
    def get_synonyms(self, cid: int) -> List[str]:
        """
        Get synonyms for a CID.
        
        Parameters:
        -----------
        cid : int
            PubChem CID
        
        Returns:
        --------
        list
            List of synonyms
        """
        syn_url = f"{self.base_url}/compound/cid/{cid}/synonyms/JSON"
        
        try:
            response = requests.get(syn_url, timeout=self.timeout)
            
            if response.status_code == 200:
                syn_data = response.json()
                if 'InformationList' in syn_data and 'Information' in syn_data['InformationList']:
                    if len(syn_data['InformationList']['Information']) > 0:
                        return syn_data['InformationList']['Information'][0].get('Synonym', [])
        except Exception as e:
            print(f"Error fetching synonyms: {str(e)}")
        
        return []
    
    def extract_chebi_id_from_synonyms(self, synonyms: List[str]) -> Optional[int]:
        """
        Extract ChEBI ID from synonym list.
        
        Parameters:
        -----------
        synonyms : list of str
            List of compound synonyms
        
        Returns:
        --------
        int or None
            ChEBI ID (without 'CHEBI:' prefix) if found, None otherwise
        
        Examples:
        ---------
        >>> client = PubChemClient()
        >>> synonyms = ['glucose', 'CHEBI:15365', 'D-glucose']
        >>> chebi_id = client.extract_chebi_id_from_synonyms(synonyms)
        >>> print(chebi_id)
        15365
        """
        # Pattern to match CHEBI:##### format
        chebi_pattern = re.compile(r'^CHEBI[:\s]+(\d+)$', re.IGNORECASE)
        
        for synonym in synonyms:
            if not isinstance(synonym, str):
                continue
            
            synonym = synonym.strip()
            match = chebi_pattern.match(synonym)
            if match:
                try:
                    return int(match.group(1))
                except ValueError:
                    continue
        
        return None
    
    def extract_chebi_ontology(self, classifications: List[Dict[str, Any]]) -> List[Any]:
        """
        Extract ChEBI ontology terms from classifications.
        
        Parameters:
        -----------
        classifications : list
            List of classification hierarchies
        
        Returns:
        --------
        list
            List of ChEBI ontology terms
        """
        chebi_ontology = []
        
        for classification in classifications:
            # Check if this is ChEBI classification
            source = classification.get('SourceName', '')
            if 'chebi' in source.lower() or 'ChEBI' in source:
                # Extract all terms from this hierarchy
                if 'Node' in classification:
                    nodes = classification['Node'] if isinstance(classification['Node'], list) else [classification['Node']]
                    for node in nodes:
                        extract_ontology_terms_from_node(node, chebi_ontology)
        
        return chebi_ontology
    
    # =========================================================================
    # Complete Compound Information
    # =========================================================================
    
    def get_compound_info_batch(
        self,
        cids: List[int]
    ) -> Dict[int, Dict[str, Any]]:
        """
        Get complete compound information for multiple CIDs.
        
        Uses efficient ChEBI ancestry-based classification when possible,
        falling back to PubChem classification if needed.
        
        Parameters:
        -----------
        cids : list of int
            List of PubChem CIDs
        
        Returns:
        --------
        dict
            Dictionary mapping CID to compound information
        """
        results = {}
        
        print(f"Processing CIDs 1-{len(cids)} of {len(cids)}...")
        
        # Get properties for all CIDs
        properties_dict = self.get_properties(cids)
        
        # Process each compound
        for cid in tqdm(cids, desc="Processing compounds", unit="compound"):
            if cid not in properties_dict:
                continue
            
            props = properties_dict[cid]
            
            # Try efficient ChEBI-based classification first
            time.sleep(self.rate_limit_delay)
            synonyms = self.get_synonyms(cid)
            chebi_id = self.extract_chebi_id_from_synonyms(synonyms)
            
            main_class = None
            subclass = None
            classifications = []
            chebi_ontology = []
            
            if chebi_id:
                # Use efficient ancestry-based classification
                try:
                    main_class, subclass = classify_by_chebi_ancestry(chebi_id)
                except Exception as e:
                    print(f"Warning: ChEBI ancestry classification failed for CID {cid} (ChEBI:{chebi_id}): {str(e)}")
                    chebi_id = None  # Force fallback
            
            # Fallback to PubChem classification if no ChEBI ID or classification failed
            if not chebi_id or main_class is None:
                time.sleep(self.rate_limit_delay)
                classifications = self.get_classification(cid)
                chebi_ontology = self.extract_chebi_ontology(classifications)
                main_class, subclass = classify_carbohydrate(chebi_ontology)
            
            is_carbohydrate = main_class is not None
            
            # Compile result
            results[cid] = {
                'pubchem_cid': cid,
                'name': props.get('IUPACName', 'Unknown'),
                'formula': props.get('MolecularFormula'),
                'molecular_weight': props.get('MolecularWeight'),
                'inchi': props.get('InChI'),
                'inchikey': props.get('InChIKey'),
                'smiles': props.get('CanonicalSMILES'),
                'chebi_id': chebi_id,
                'classifications': classifications,
                'chebi_ontology': chebi_ontology,
                'is_carbohydrate': is_carbohydrate,
                'carbohydrate_main_class': main_class,
                'carbohydrate_subclass': subclass
            }
        
        return results


# =============================================================================
# High-level API Functions
# =============================================================================

# Global client instance
_default_client = PubChemClient()


def get_compound_info_pubchem(
    identifier: Union[str, List[str]],
    identifier_type: str = 'auto'
) -> Union[Optional[Dict[str, Any]], List[Optional[Dict[str, Any]]]]:
    """
    Fetch compound information from PubChem using InChIKey or SMILES.
    
    Supports both single identifier and batch queries (list of identifiers).
    Always uses POST API for efficiency.
    
    Parameters:
    -----------
    identifier : str or list of str
        Single identifier or list of identifiers (InChIKey or SMILES strings)
        When providing a list, all identifiers must be of the same type
    identifier_type : str, optional
        Type of identifier: 'inchikey', 'smiles', or 'auto' (default: 'auto')
        'auto' will detect based on format of first identifier
    
    Returns:
    --------
    dict or list of dict or None
        For single identifier: Dictionary or None
        For list of identifiers: List of dictionaries (None for failed queries)
        
        Dictionary containing compound information:
        - pubchem_cid: PubChem Compound ID
        - name: IUPAC name
        - formula: Molecular formula
        - molecular_weight: Molecular weight
        - inchi: InChI string
        - inchikey: InChIKey
        - smiles: Canonical SMILES
        - chebi_id: ChEBI ID (int) if found in synonyms, None otherwise
        - classifications: List of chemical classification hierarchies (empty if ChEBI method used)
        - chebi_ontology: ChEBI ontology terms extracted from classifications (empty if ChEBI method used)
        - is_carbohydrate: Boolean flag (True if compound belongs to CHEBI:78616)
        - carbohydrate_main_class: Main classification category (5 categories or None)
        - carbohydrate_subclass: Subclass name or None
    
    Example:
    --------
    >>> # Single query
    >>> result = get_compound_info_pubchem("RBNPOMFGQQGHHO-UHFFFAOYSA-N")
    >>> print(f"Compound: {result['name']}, Is carbohydrate: {result['is_carbohydrate']}")
    >>> 
    >>> # Batch query
    >>> identifiers = ["RBNPOMFGQQGHHO-UHFFFAOYSA-N", "WQZGKKKJIJFFOK-GASJEMHNSA-N"]
    >>> results = get_compound_info_pubchem(identifiers)
    >>> for result in results:
    >>>     if result:
    >>>         print(f"Compound: {result['name']}")
    """
    client = _default_client
    
    # Determine if input is single or list
    is_single = isinstance(identifier, str)
    
    # Convert single identifier to list for unified processing
    identifiers = [identifier] if is_single else identifier
    
    if not identifiers:
        return None if is_single else []
    
    # Auto-detect or validate identifier type
    first_id = identifiers[0]
    if identifier_type == 'auto':
        if '-' in first_id and len(first_id.split('-')) == 3:
            identifier_type = 'inchikey'
        else:
            identifier_type = 'smiles'
    
    # Validate all identifiers are of the same type
    for idx, id_str in enumerate(identifiers):
        if identifier_type == 'inchikey':
            if not ('-' in id_str and len(id_str.split('-')) == 3):
                raise ValueError(
                    f"Identifier at index {idx} ('{id_str}') is not a valid InChIKey. "
                    f"All identifiers must be of the same type."
                )
        else:  # smiles
            if '-' in id_str and len(id_str.split('-')) == 3:
                raise ValueError(
                    f"Identifier at index {idx} ('{id_str}') appears to be an InChIKey "
                    f"but SMILES was expected. All identifiers must be of the same type."
                )
    
    # Process using POST API
    if len(identifiers) == 1:
        print(f"Processing single {identifier_type} identifier...")
    else:
        print(f"Processing batch of {len(identifiers)} {identifier_type} identifiers...")
    
    # Step 1: Resolve identifiers to CIDs
    identifier_to_cid = client.resolve_identifiers_to_cids(identifiers, identifier_type)
    
    # Step 2: Get all valid CIDs
    valid_cids = [cid for cid in identifier_to_cid.values() if cid is not None]
    
    if not valid_cids:
        print("No valid CIDs found")
        if is_single:
            return None
        else:
            return [None for _ in identifiers]
    
    print(f"\nFound {len(valid_cids)} valid CIDs, fetching properties...")
    
    # Step 3: Batch fetch complete information
    results_dict = client.get_compound_info_batch(valid_cids)
    
    # Step 4: Create result list in original order
    results = []
    for identifier in identifiers:
        cid = identifier_to_cid.get(identifier)
        if cid and cid in results_dict:
            results.append(results_dict[cid])
        else:
            results.append(None)
    
    # Return single result or list based on input type
    return results[0] if is_single else results
