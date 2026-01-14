"""
ChEBI API Client

This module provides functions for interacting with the ChEBI (Chemical Entities of
Biological Interest) API to retrieve ontology information about chemical compounds.
"""

import requests
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from . import config

# Cache for ChEBI API responses to avoid repeated calls
_chebi_children_cache: Dict[int, List[Dict[str, Any]]] = {}


class ChEBIClient:
    """
    Client for interacting with ChEBI API.
    
    This class provides methods to query the ChEBI ontology database
    and retrieve hierarchical relationships between chemical entities.
    Supports persistent caching to reduce API calls.
    """
    
    def __init__(self, base_url: Optional[str] = None, timeout: Optional[int] = None, 
                 cache_dir: Optional[Path] = None, use_cache: bool = True):
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
        
        self.cache_file = self.cache_dir / 'chebi_children_cache.json'
        
        # Use global cache for backward compatibility
        self.cache = _chebi_children_cache
        
        # Track if cache has been modified
        self._cache_dirty = False
        self._cache_save_counter = 0
        self._save_batch_size = 10  # Save after every N new entries
        
        # Load persistent cache if enabled
        if self.use_cache:
            self._load_cache()
    
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
        if chebi_id in self.cache:
            return self.cache[chebi_id]
        
        try:
            url = f'{self.base_url}/ontology/children/{chebi_id}/'
            response = requests.get(url, headers={'accept': '*/*'}, timeout=self.timeout)
            
            if response.status_code != 200:
                print(f"Warning: ChEBI API failed for {chebi_id}: HTTP {response.status_code}")
                self.cache[chebi_id] = []
                return []
            
            data = response.json()
            incoming = data.get('ontology_relations', {}).get('incoming_relations', [])
            
            # Filter for "is a" relations (direct children)
            children = [rel for rel in incoming if rel.get('relation_type') == 'is a']
            
            self.cache[chebi_id] = children
            self._cache_dirty = True
            self._cache_save_counter += 1
            
            # Save cache periodically (batch saves)
            if self.use_cache and self._cache_save_counter >= self._save_batch_size:
                self._save_cache()
                self._cache_save_counter = 0
            
            return children
            
        except Exception as e:
            print(f"Warning: Error fetching ChEBI children for {chebi_id}: {str(e)}")
            self.cache[chebi_id] = []
            self._cache_dirty = True
            return []
    
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
    
    def clear_cache(self) -> None:
        """Clear the ChEBI children cache."""
        self.cache.clear()
    
    def _load_cache(self) -> None:
        """
        Load cache from disk if it exists.
        
        The cache is stored as JSON with ChEBI IDs as string keys
        (JSON doesn't support integer keys).
        """
        if not self.cache_file.exists():
            return
        
        try:
            with open(self.cache_file, 'r', encoding='utf-8') as f:
                cache_data = json.load(f)
            
            # Convert string keys back to integers
            for key_str, value in cache_data.items():
                self.cache[int(key_str)] = value
            
            print(f"Loaded ChEBI cache with {len(self.cache)} entries from {self.cache_file}")
        except Exception as e:
            print(f"Warning: Failed to load ChEBI cache: {str(e)}")
    
    def _save_cache(self) -> None:
        """
        Save cache to disk in JSON format.
        
        Creates the cache directory if it doesn't exist.
        Converts integer keys to strings for JSON compatibility.
        """
        if not self.use_cache or not self._cache_dirty:
            return
        
        try:
            # Create cache directory if it doesn't exist
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            
            # Convert integer keys to strings for JSON
            cache_data = {str(key): value for key, value in self.cache.items()}
            
            # Write to temporary file first, then rename for atomicity
            temp_file = self.cache_file.with_suffix('.tmp')
            with open(temp_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, indent=2, ensure_ascii=False)
            
            # Atomic rename
            temp_file.replace(self.cache_file)
            
            self._cache_dirty = False
            print(f"Saved ChEBI cache with {len(self.cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to save ChEBI cache: {str(e)}")
    
    def save_cache(self) -> None:
        """
        Manually save the cache to disk.
        
        This is called automatically in batches and on cleanup,
        but can be called manually to ensure data is persisted.
        """
        self._save_cache()
    
    def __del__(self):
        """
        Destructor to save cache when object is garbage collected.
        
        This ensures that any pending cache updates are saved even
        if save_cache() wasn't called explicitly.
        """
        if self.use_cache and self._cache_dirty:
            self._save_cache()


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
