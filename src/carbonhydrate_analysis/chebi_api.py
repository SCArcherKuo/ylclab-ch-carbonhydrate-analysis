"""
ChEBI API Client

This module provides functions for interacting with the ChEBI (Chemical Entities of
Biological Interest) API to retrieve ontology information about chemical compounds.
"""

import requests
from pathlib import Path
from typing import List, Dict, Any, Optional
from loguru import logger
from . import config
from .cache_manager import CacheManager, PersistentCache
from .retry import retry_with_backoff


class ChEBIClient:
    """
    Client for interacting with ChEBI API.
    
    This class provides methods to query the ChEBI ontology database
    and retrieve hierarchical relationships between chemical entities.
    Supports persistent caching to reduce API calls.
    """
    
    def __init__(
        self,
        base_url: Optional[str] = None,
        timeout: Optional[int] = None,
        cache_dir: Optional[Path] = None,
        use_cache: bool = True
    ):
        """
        Initialize ChEBI API client.
        
        Parameters:
        -----------
        base_url : str, optional
            Base URL for ChEBI API (defaults to config.CHEBI_API_BASE_URL)
        timeout : int, optional
            Request timeout in seconds (defaults to config.API_TIMEOUT)
        cache_dir : Path, optional
            Directory for cache storage (defaults to data/cache/)
        use_cache : bool, optional
            Whether to use persistent caching (default: True)
        """
        self.base_url = base_url or config.CHEBI_API_BASE_URL
        self.timeout = timeout or config.API_TIMEOUT
        self.use_cache = use_cache
        
        # Set up cache directory
        if cache_dir:
            self.cache_dir = cache_dir
        else:
            # Default to data/cache relative to project root
            project_root = Path(__file__).parent.parent.parent
            self.cache_dir = project_root / 'data' / 'cache'
        
        # Initialize cache manager
        if self.use_cache:
            self.cache_manager = CacheManager(self.cache_dir)
            self.children_cache = self.cache_manager.register_cache(
                'chebi_children_cache',
                cache_file=self.cache_dir / 'chebi_children_cache.json'
            )
            self.parents_cache = self.cache_manager.register_cache(
                'chebi_parents_cache',
                cache_file=self.cache_dir / 'chebi_parents_cache.json'
            )
            self.ancestors_cache = self.cache_manager.register_cache(
                'chebi_ancestors_cache',
                cache_file=self.cache_dir / 'chebi_ancestors_cache.json'
            )
        else:
            self.cache_manager = None
            self.children_cache = None
            self.parents_cache = None
            self.ancestors_cache = None
    
    @retry_with_backoff(
        max_retries=3,
        base_delay=1.0,
        exceptions=(requests.exceptions.Timeout, requests.exceptions.RequestException)
    )
    def _fetch_children(self, chebi_id: int) -> List[Dict[str, Any]]:
        """
        Fetch children from ChEBI API with retry logic.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID (without 'CHEBI:' prefix)
        
        Returns:
        --------
        list of dict
            List of direct children with 'is a' relationship
        """
        url = f'{self.base_url}/ontology/children/{chebi_id}/'
        response = requests.get(url, headers={'accept': '*/*'}, timeout=self.timeout)
        
        if response.status_code != 200:
            logger.warning(f"ChEBI API failed for {chebi_id}: HTTP {response.status_code}")
            return []
        
        data = response.json()
        incoming = data.get('ontology_relations', {}).get('incoming_relations', [])
        
        # Filter for "is a" relations (direct children)
        children = [rel for rel in incoming if rel.get('relation_type') == 'is a']
        return children
    
    def get_children(self, chebi_id: int) -> List[Dict[str, Any]]:
        """
        Get direct children of a ChEBI ID from ChEBI backend API.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID (without 'CHEBI:' prefix)
        
        Returns:
        --------
        list of dict
            List of direct children with 'is a' relationship.
            Each dict contains: 'init_id', 'init_name', 'relation_type'
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> children = client.get_children(16646)  # carbohydrate
        >>> print(len(children))
        """
        # Check cache first
        if self.use_cache and chebi_id in self.children_cache:
            logger.debug(f"Cache hit for ChEBI children: {chebi_id}")
            return self.children_cache[chebi_id]
        
        logger.debug(f"Fetching ChEBI children for {chebi_id}")
        
        try:
            children = self._fetch_children(chebi_id)
            logger.debug(f"Found {len(children)} children for ChEBI {chebi_id}")
            
            # Cache the result
            if self.use_cache and children:  # Only cache non-empty results
                self.children_cache[chebi_id] = children
            
            return children
            
        except Exception as e:
            logger.error(f"Error fetching ChEBI children for {chebi_id}: {str(e)}")
            return []
    
    @retry_with_backoff(
        max_retries=3,
        base_delay=1.0,
        exceptions=(requests.exceptions.Timeout, requests.exceptions.RequestException)
    )
    def _fetch_parents(self, chebi_id: int) -> List[Dict[str, Any]]:
        """
        Fetch parents from ChEBI API with retry logic.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID (without 'CHEBI:' prefix)
        
        Returns:
        --------
        list of dict
            List of direct parents with 'is a' relationship
        """
        url = f'{self.base_url}/ontology/parents/{chebi_id}/'
        response = requests.get(url, headers={'accept': '*/*'}, timeout=self.timeout)
        
        if response.status_code != 200:
            logger.warning(f"ChEBI API failed for {chebi_id}: HTTP {response.status_code}")
            return []
        
        data = response.json()
        outgoing = data.get('ontology_relations', {}).get('outgoing_relations', [])
        
        # Filter for "is a" relations (direct parents)
        parents = [rel for rel in outgoing if rel.get('relation_type') == 'is a']
        return parents
    
    def get_parents(self, chebi_id: int) -> List[Dict[str, Any]]:
        """
        Get direct parents of a ChEBI ID from ChEBI backend API.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID (without 'CHEBI:' prefix)
        
        Returns:
        --------
        list of dict
            List of direct parents with 'is a' relationship.
            Each dict contains: 'init_id', 'init_name', 'relation_type', 'final_id', 'final_name'
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> parents = client.get_parents(4167)  # D-glucopyranose
        >>> print(parents)
        """
        # Check cache first
        if self.use_cache and chebi_id in self.parents_cache:
            logger.debug(f"Cache hit for ChEBI parents: {chebi_id}")
            return self.parents_cache[chebi_id]
        
        logger.debug(f"Fetching ChEBI parents for {chebi_id}")
        
        try:
            parents = self._fetch_parents(chebi_id)
            logger.debug(f"Found {len(parents)} parents for ChEBI {chebi_id}")
            
            # Cache the result
            if self.use_cache and parents:  # Only cache non-empty results
                self.parents_cache[chebi_id] = parents
            
            return parents
            
        except Exception as e:
            logger.error(f"Error fetching ChEBI parents for {chebi_id}: {str(e)}")
            return []
    
    def get_all_ancestors(self, chebi_id: int, max_depth: int = 20) -> List[int]:
        """
        Get all ancestor ChEBI IDs by recursively traversing parent relationships.
        
        Parameters:
        -----------
        chebi_id : int
            Starting ChEBI ID (without 'CHEBI:' prefix)
        max_depth : int, optional
            Maximum recursion depth to prevent infinite loops (default: 20)
        
        Returns:
        --------
        list of int
            List of all ancestor ChEBI IDs (includes the starting ID)
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> ancestors = client.get_all_ancestors(4167)  # D-glucopyranose
        >>> print(16646 in ancestors)  # Check if carbohydrate ancestor
        True
        """
        # Check cache first
        if self.use_cache and chebi_id in self.ancestors_cache:
            logger.debug(f"Cache hit for ChEBI ancestors: {chebi_id}")
            return self.ancestors_cache[chebi_id]
        
        logger.debug(f"Computing all ancestors for ChEBI {chebi_id}")
        ancestors = set([chebi_id])  # Include self
        visited = set()
        queue = [chebi_id]
        depth = 0
        
        try:
            while queue and depth < max_depth:
                current_batch = queue[:]
                queue = []
                
                for current_id in current_batch:
                    if current_id in visited:
                        continue
                    visited.add(current_id)
                    
                    parents = self.get_parents(current_id)
                    for parent_rel in parents:
                        parent_id = parent_rel.get('final_id')
                        if parent_id and parent_id not in ancestors:
                            ancestors.add(parent_id)
                            queue.append(parent_id)
                
                depth += 1
            
            ancestor_list = sorted(list(ancestors))
            logger.debug(f"Found {len(ancestor_list)} ancestors for ChEBI {chebi_id}")
            
            # Cache the result
            if self.use_cache:
                self.ancestors_cache[chebi_id] = ancestor_list
            
            return ancestor_list
            
        except Exception as e:
            logger.exception(f"Error computing ancestors for ChEBI {chebi_id}: {str(e)}")
            # Return partial results if available
            ancestor_list = sorted(list(ancestors))
            return ancestor_list
    
    def get_main_groups(self, chebi_id: int) -> List[str]:
        """
        Get main groups (children with >1 children) of a ChEBI ID.
        
        Main groups are defined as direct children that themselves have
        more than one child, indicating significant branching in the ontology.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID to get main groups for
        
        Returns:
        --------
        list of str
            List of main group names (children that have >1 children themselves)
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> main_groups = client.get_main_groups(16646)  # carbohydrate
        >>> print(main_groups)
        """
        children = self.get_children(chebi_id)
        main_groups = []
        
        for child in children:
            child_id = child.get('init_id')
            child_name = child.get('init_name', '')
            
            if child_id:
                # Check if this child has more than 1 children
                grandchildren = self.get_children(child_id)
                if len(grandchildren) > 1:
                    main_groups.append(child_name)
        
        return main_groups
    
    def save_cache(self) -> None:
        """Manually save all caches to disk."""
        if self.use_cache and self.cache_manager:
            self.cache_manager.save_all()
    
    def clear_cache(self) -> None:
        """Clear all caches."""
        if self.use_cache and self.cache_manager:
            self.cache_manager.clear_all()
    
    def __del__(self):
        """Destructor to save caches when object is garbage collected."""
        if self.use_cache and self.cache_manager:
            self.cache_manager.save_all()


# =============================================================================
# Module-level convenience functions for backward compatibility
# =============================================================================

# Global client instance
_default_client = ChEBIClient()


def get_chebi_children(chebi_id: int) -> List[Dict[str, Any]]:
    """
    Get direct children of a ChEBI ID from ChEBI backend API.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of dict
        List of direct children with 'is a' relationship.
    """
    return _default_client.get_children(chebi_id)


def get_main_groups(chebi_id: int) -> List[str]:
    """
    Get main groups (children with >1 children) of a ChEBI ID.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID to get main groups for
    
    Returns:
    --------
    list of str
        List of main group names
    """
    return _default_client.get_main_groups(chebi_id)


def get_chebi_parents(chebi_id: int) -> List[Dict[str, Any]]:
    """
    Get direct parents of a ChEBI ID from ChEBI backend API.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of dict
        List of direct parents with 'is a' relationship.
    """
    return _default_client.get_parents(chebi_id)


def get_all_ancestors(chebi_id: int) -> List[int]:
    """
    Get all ancestor ChEBI IDs by recursively traversing parent relationships.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        Starting ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of int
        List of all ancestor ChEBI IDs (includes the starting ID)
    """
    return _default_client.get_all_ancestors(chebi_id)
