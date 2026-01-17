"""
PubChem API Client (Refactored)

This module provides functions for interacting with the PubChem API to retrieve
compound information, properties, classifications, and synonyms.
"""

import requests
import json
import time
import re
from pathlib import Path
from typing import List, Dict, Any, Union, Optional
from tqdm import tqdm
from loguru import logger
from . import config
from .classification import classify_carbohydrate, classify_by_chebi_ancestry
from .utils import extract_ontology_terms_from_node
from .retry import RetryManager
from .rate_limiter import RateLimiter
from .error_tracker import ServerErrorTracker, FailedIdentifierTracker


class PubChemClient:
    """
    Client for interacting with PubChem API.
    
    This class provides methods to query PubChem database for compound
    information using various identifiers (InChIKey, SMILES, CID).
    """
    
    def __init__(
        self,
        base_url: Optional[str] = None,
        timeout: Optional[int] = None,
        rate_limit_delay: Optional[float] = None
    ):
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
        self.chunk_size = 512
        self.max_synonyms = config.MAX_SYNONYMS
        
        # Initialize utility components
        self.retry_manager = RetryManager(
            max_retries=config.MAX_RETRIES,
            base_delay=config.RETRY_DELAY,
            backoff_factor=config.RETRY_BACKOFF
        )
        self.rate_limiter = RateLimiter(delay=rate_limit_delay or config.RATE_LIMIT_DELAY)
        self.error_tracker = ServerErrorTracker(
            error_window=config.SERVER_ERROR_WINDOW,
            error_threshold=config.SERVER_ERROR_THRESHOLD,
            sleep_time=config.SERVER_ERROR_SLEEP_TIME
        )
        
        # Set up cache directory for failed identifiers
        project_root = Path(__file__).parent.parent.parent
        cache_dir = project_root / 'data' / 'cache'
        self.failed_tracker = FailedIdentifierTracker(cache_dir)
    
    # =========================================================================
    # Identifier Resolution
    # =========================================================================
    
    def resolve_identifier_to_cid(self, identifier: str, identifier_type: str) -> Optional[int]:
        """
        Resolve a single identifier to PubChem CID with retry logic.
        
        Parameters:
        -----------
        identifier : str
            Identifier to resolve (InChIKey or SMILES)
        identifier_type : str
            Type of identifier: 'inchikey' or 'smiles'
        
        Returns:
        --------
        int or None
            PubChem CID if found, None otherwise
        """
        logger.debug(f"Resolving {identifier_type}: {identifier[:30]}...")
        
        def _resolve():
            url = f'{self.base_url}/compound/{identifier_type}/{identifier}/cids/JSON'
            self.rate_limiter.wait()
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                    cids = data['IdentifierList']['CID']
                    if cids:
                        cid = cids[0]
                        logger.debug(f"Resolved {identifier[:30]}... to CID {cid}")
                        return cid
                logger.warning(f"No CID found for {identifier[:30]}...")
                return None
            elif response.status_code == 429:
                logger.error(f"Rate limit exceeded for {identifier[:30]}...")
                return None
            elif response.status_code >= 500:
                self.error_tracker.record_error()
                raise requests.exceptions.RequestException(f"Server error: HTTP {response.status_code}")
            else:
                logger.warning(f"Identifier not found: {identifier[:30]}... (HTTP {response.status_code})")
                return None
        
        try:
            return self.retry_manager.execute_with_retry(
                _resolve,
                exceptions=(requests.exceptions.Timeout, requests.exceptions.RequestException)
            )
        except Exception as e:
            logger.error(f"Failed to resolve {identifier_type} {identifier[:30]}...: {str(e)}")
            self.failed_tracker.add_failed_identifier(identifier_type, identifier)
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
        
        # Save any failed identifiers
        self.failed_tracker.save_failed_identifiers(identifier_type)
        
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
        logger.info(f"Fetching properties for {len(cids)} CIDs in chunks of {self.chunk_size}")
        properties = {}
        
        for i in range(0, len(cids), self.chunk_size):
            chunk_cids = cids[i:i + self.chunk_size]
            chunk_num = i // self.chunk_size + 1
            total_chunks = (len(cids) + self.chunk_size - 1) // self.chunk_size
            logger.debug(f"Processing chunk {chunk_num}/{total_chunks} ({len(chunk_cids)} CIDs)")
            
            chunk_properties = self._fetch_properties_chunk(chunk_cids, chunk_num)
            properties.update(chunk_properties)
        
        logger.info(f"Fetched properties for {len(properties)}/{len(cids)} CIDs")
        return properties
    
    def _fetch_properties_chunk(self, cids: List[int], chunk_num: int) -> Dict[int, Dict[str, Any]]:
        """Fetch properties for a single chunk of CIDs."""
        cid_list = ','.join(map(str, cids))
        
        def _fetch():
            self.rate_limiter.wait()
            props_url = f"{self.base_url}/compound/cid/property/MolecularFormula,MolecularWeight,InChI,InChIKey,SMILES,IUPACName/JSON"
            response = requests.post(
                props_url,
                data=f"cid={cid_list}",
                headers={'Content-Type': 'application/x-www-form-urlencoded'},
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                props_data = response.json()
                properties = {}
                if 'PropertyTable' in props_data and 'Properties' in props_data['PropertyTable']:
                    for props in props_data['PropertyTable']['Properties']:
                        cid = props.get('CID')
                        if cid:
                            properties[cid] = props
                logger.debug(f"Successfully fetched properties for chunk {chunk_num}")
                return properties
            elif response.status_code == 429:
                logger.error(f"Rate limit exceeded for properties chunk {chunk_num}")
                return {}
            elif response.status_code >= 500:
                self.error_tracker.record_error()
                raise requests.exceptions.RequestException(f"Server error: HTTP {response.status_code}")
            else:
                logger.warning(f"Properties request failed for chunk {chunk_num}: HTTP {response.status_code}")
                return {}
        
        try:
            return self.retry_manager.execute_with_retry(
                _fetch,
                exceptions=(requests.exceptions.Timeout, requests.exceptions.RequestException)
            )
        except Exception as e:
            logger.error(f"Failed to fetch properties for chunk {chunk_num}: {str(e)}")
            for cid in cids:
                self.failed_tracker.add_failed_cid(cid)
            return {}
    
    # =========================================================================
    # Classification and Synonyms
    # =========================================================================
    
    def get_classification(self, cid: int) -> List[Dict[str, Any]]:
        """
        Get classification hierarchies for a CID with retry logic.
        
        Parameters:
        -----------
        cid : int
            PubChem CID
        
        Returns:
        --------
        list
            List of classification hierarchies
        """
        logger.debug(f"Fetching classification for CID {cid}")
        
        def _fetch():
            self.rate_limiter.wait()
            class_url = f"{self.base_url}/compound/cid/{cid}/classification/JSON"
            response = requests.get(class_url, timeout=self.timeout)
            
            if response.status_code == 200:
                class_data = response.json()
                if 'Hierarchies' in class_data and 'Hierarchy' in class_data['Hierarchies']:
                    hierarchies = [h for h in class_data['Hierarchies']['Hierarchy'] if 'Node' in h]
                    logger.debug(f"Found {len(hierarchies)} classification hierarchies for CID {cid}")
                    return hierarchies
                return []
            elif response.status_code == 429:
                logger.error(f"Rate limit exceeded for classification CID {cid}")
                return []
            elif response.status_code >= 500:
                self.error_tracker.record_error()
                raise requests.exceptions.RequestException(f"Server error: HTTP {response.status_code}")
            else:
                logger.warning(f"No classification data for CID {cid}: HTTP {response.status_code}")
                return []
        
        try:
            return self.retry_manager.execute_with_retry(
                _fetch,
                exceptions=(requests.exceptions.Timeout, requests.exceptions.RequestException)
            )
        except Exception as e:
            logger.error(f"Failed to fetch classification for CID {cid}: {str(e)}")
            self.failed_tracker.add_failed_cid(cid)
            return []
    
    def get_synonyms(self, cid: int) -> List[str]:
        """
        Get synonyms for a CID with retry logic.
        
        Parameters:
        -----------
        cid : int
            PubChem CID
        
        Returns:
        --------
        list
            List of synonyms
        """
        logger.debug(f"Fetching synonyms for CID {cid}")
        
        def _fetch():
            self.rate_limiter.wait()
            syn_url = f"{self.base_url}/compound/cid/{cid}/synonyms/JSON"
            response = requests.get(syn_url, timeout=self.timeout)
            
            if response.status_code == 200:
                syn_data = response.json()
                if 'InformationList' in syn_data and 'Information' in syn_data['InformationList']:
                    if len(syn_data['InformationList']['Information']) > 0:
                        synonyms = syn_data['InformationList']['Information'][0].get('Synonym', [])
                        logger.debug(f"Found {len(synonyms)} synonyms for CID {cid}")
                        return synonyms
                return []
            elif response.status_code == 429:
                logger.error(f"Rate limit exceeded for synonyms CID {cid}")
                return []
            elif response.status_code >= 500:
                self.error_tracker.record_error()
                raise requests.exceptions.RequestException(f"Server error: HTTP {response.status_code}")
            else:
                logger.warning(f"No synonyms found for CID {cid}: HTTP {response.status_code}")
                return []
        
        try:
            return self.retry_manager.execute_with_retry(
                _fetch,
                exceptions=(requests.exceptions.Timeout, requests.exceptions.RequestException)
            )
        except Exception as e:
            logger.error(f"Failed to fetch synonyms for CID {cid}: {str(e)}")
            self.failed_tracker.add_failed_cid(cid)
            return []
    
    # =========================================================================
    # Data Extraction Helpers
    # =========================================================================
    
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
        """
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
            source = classification.get('SourceName', '')
            if 'chebi' in source.lower() or 'ChEBI' in source:
                if 'Node' in classification:
                    nodes = classification['Node'] if isinstance(classification['Node'], list) else [classification['Node']]
                    for node in nodes:
                        extract_ontology_terms_from_node(node, chebi_ontology)
        
        return chebi_ontology
    
    # =========================================================================
    # Complete Compound Information
    # =========================================================================
    
    def get_compound_info_batch(self, cids: List[int]) -> Dict[int, Dict[str, Any]]:
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
        logger.info(f"Starting batch processing of {len(cids)} compounds")
        results = {}
        
        print(f"Processing CIDs 1-{len(cids)} of {len(cids)}...")
        
        # Get properties for all CIDs
        properties_dict = self.get_properties(cids)
        logger.info(f"Retrieved properties for {len(properties_dict)} compounds")
        
        # Process each compound
        for idx, cid in enumerate(tqdm(cids, desc="Processing compounds", unit="compound"), 1):
            try:
                result = self._process_single_compound(cid, idx, len(cids), properties_dict)
                if result:
                    results[cid] = result
            except KeyboardInterrupt:
                logger.warning(f"Keyboard interrupt received at compound {idx}/{len(cids)} (CID {cid})")
                raise
            except Exception as e:
                logger.exception(f"Unexpected error processing CID {cid} at position {idx}/{len(cids)}: {str(e)}")
                continue
        
        logger.info(f"Batch processing complete: {len(results)}/{len(cids)} compounds processed successfully")
        
        # Save any failed CIDs
        self.failed_tracker.save_failed_cids()
        
        return results
    
    def _process_single_compound(
        self,
        cid: int,
        idx: int,
        total: int,
        properties_dict: Dict[int, Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """Process a single compound and return its information."""
        logger.debug(f"Processing compound {idx}/{total}: CID {cid}")
        
        if cid not in properties_dict:
            logger.warning(f"CID {cid} not in properties dict, skipping")
            return None
        
        props = properties_dict[cid]
        logger.debug(f"CID {cid}: {props.get('IUPACName', 'Unknown')[:50]}...")
        
        # Try efficient ChEBI-based classification first
        logger.debug(f"CID {cid}: Fetching synonyms")
        synonyms = self.get_synonyms(cid)
        chebi_id = self.extract_chebi_id_from_synonyms(synonyms)
        
        if chebi_id:
            logger.debug(f"CID {cid}: Found ChEBI ID {chebi_id}")
        
        main_class = None
        subclass = None
        classifications = []
        chebi_ontology = []
        
        if chebi_id:
            # Use efficient ancestry-based classification
            try:
                logger.debug(f"CID {cid}: Classifying via ChEBI ancestry")
                main_class, subclass = classify_by_chebi_ancestry(chebi_id)
                logger.debug(f"CID {cid}: ChEBI classification result: {main_class}, {subclass}")
            except Exception as e:
                logger.error(f"ChEBI ancestry classification failed for CID {cid} (ChEBI:{chebi_id}): {str(e)}")
                chebi_id = None  # Force fallback on error
        
        # Fallback to PubChem classification ONLY if no ChEBI ID found
        if not chebi_id:
            logger.debug(f"CID {cid}: No ChEBI ID found, falling back to PubChem classification")
            classifications = self.get_classification(cid)
            chebi_ontology = self.extract_chebi_ontology(classifications)
            main_class, subclass = classify_carbohydrate(chebi_ontology)
            logger.debug(f"CID {cid}: PubChem classification result: {main_class}, {subclass}")
        
        is_carbohydrate = main_class is not None
        
        # Compile result
        return {
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
        return None if is_single else [None for _ in identifiers]
    
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
