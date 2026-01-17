"""
Cache Manager

This module provides a unified caching system with LRU in-memory caching
and persistent JSON file storage.
"""

import json
from pathlib import Path
from typing import Any, Dict, Optional, TypeVar, Generic
from collections import OrderedDict
from loguru import logger

T = TypeVar('T')


class LRUCache(Generic[T]):
    """
    Least Recently Used (LRU) cache with maximum size limit.
    
    Automatically evicts least recently used items when cache is full.
    This prevents unbounded memory growth during long-running processes.
    """
    
    def __init__(self, max_size: int = 200):
        """
        Initialize LRU cache.
        
        Parameters:
        -----------
        max_size : int
            Maximum number of items to keep in cache (default: 200)
        """
        self.cache: OrderedDict[Any, T] = OrderedDict()
        self.max_size = max_size
    
    def get(self, key: Any, default: Optional[T] = None) -> Optional[T]:
        """
        Get value from cache, moving it to end (most recently used).
        
        Parameters:
        -----------
        key : any
            Cache key
        default : any, optional
            Default value if key not found
        
        Returns:
        --------
        any
            Cached value or default
        """
        if key in self.cache:
            # Move to end (most recent)
            self.cache.move_to_end(key)
            return self.cache[key]
        return default
    
    def __getitem__(self, key: Any) -> T:
        """Dictionary-style access."""
        value = self.get(key)
        if value is None:
            raise KeyError(key)
        return value
    
    def __setitem__(self, key: Any, value: T) -> None:
        """Set item in cache."""
        self.set(key, value)
    
    def set(self, key: Any, value: T) -> None:
        """
        Set value in cache.
        
        Parameters:
        -----------
        key : any
            Cache key
        value : any
            Value to cache
        """
        if key in self.cache:
            # Update existing key
            self.cache.move_to_end(key)
        else:
            # Add new key
            if len(self.cache) >= self.max_size:
                # Remove least recently used
                self.cache.popitem(last=False)
        
        self.cache[key] = value
    
    def __contains__(self, key: Any) -> bool:
        """Check if key exists in cache."""
        return key in self.cache
    
    def clear(self) -> None:
        """Clear all cached items."""
        self.cache.clear()
    
    def size(self) -> int:
        """Get current cache size."""
        return len(self.cache)
    
    def keys(self):
        """Get cache keys."""
        return self.cache.keys()
    
    def values(self):
        """Get cache values."""
        return self.cache.values()
    
    def items(self):
        """Get cache items."""
        return self.cache.items()


class PersistentCache(Generic[T]):
    """
    Cache with persistent JSON file storage.
    
    Combines in-memory LRU cache with disk persistence for long-term storage.
    Supports automatic saving at configurable intervals.
    """
    
    def __init__(
        self,
        cache_file: Path,
        max_memory_size: int = 200,
        save_batch_size: int = 10,
        auto_save: bool = True
    ):
        """
        Initialize persistent cache.
        
        Parameters:
        -----------
        cache_file : Path
            Path to JSON cache file
        max_memory_size : int
            Maximum items in memory cache
        save_batch_size : int
            Save to disk after this many new entries
        auto_save : bool
            Whether to automatically save periodically
        """
        self.cache_file = cache_file
        self.memory_cache: LRUCache[T] = LRUCache(max_size=max_memory_size)
        self.save_batch_size = save_batch_size
        self.auto_save = auto_save
        self.save_counter = 0
        self.is_dirty = False
        
        # Ensure cache directory exists
        self.cache_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Load from disk
        self.load()
    
    def get(self, key: Any, default: Optional[T] = None) -> Optional[T]:
        """
        Get value from cache.
        
        Parameters:
        -----------
        key : any
            Cache key
        default : any, optional
            Default value if key not found
        
        Returns:
        --------
        any
            Cached value or default
        """
        return self.memory_cache.get(key, default)
    
    def set(self, key: Any, value: T) -> None:
        """
        Set value in cache.
        
        Parameters:
        -----------
        key : any
            Cache key
        value : any
            Value to cache
        """
        was_new = key not in self.memory_cache
        self.memory_cache.set(key, value)
        
        if was_new:
            self.is_dirty = True
            self.save_counter += 1
            
            # Auto-save if threshold reached
            if self.auto_save and self.save_counter >= self.save_batch_size:
                self.save()
                self.save_counter = 0
    
    def __getitem__(self, key: Any) -> T:
        """Dictionary-style access."""
        return self.memory_cache[key]
    
    def __setitem__(self, key: Any, value: T) -> None:
        """Set item in cache."""
        self.set(key, value)
    
    def __contains__(self, key: Any) -> bool:
        """Check if key exists in cache."""
        return key in self.memory_cache
    
    def load(self) -> None:
        """Load cache from disk."""
        if not self.cache_file.exists():
            logger.debug(f"Cache file not found: {self.cache_file}")
            return
        
        try:
            with open(self.cache_file, 'r') as f:
                data = json.load(f)
            
            # Convert string keys back to appropriate types if needed
            for key, value in data.items():
                # Try to convert key to int if possible
                try:
                    key = int(key)
                except (ValueError, TypeError):
                    pass
                self.memory_cache.set(key, value)
            
            logger.info(f"Loaded {len(data)} entries from {self.cache_file.name}")
        except Exception as e:
            logger.error(f"Failed to load cache from {self.cache_file}: {e}")
    
    def save(self) -> None:
        """Save cache to disk."""
        if not self.is_dirty:
            return
        
        try:
            # Convert cache to dict for JSON serialization
            data = dict(self.memory_cache.items())
            
            with open(self.cache_file, 'w') as f:
                json.dump(data, f, indent=2)
            
            self.is_dirty = False
            logger.debug(f"Saved {len(data)} entries to {self.cache_file.name}")
        except Exception as e:
            logger.error(f"Failed to save cache to {self.cache_file}: {e}")
    
    def clear(self) -> None:
        """Clear cache in memory and on disk."""
        self.memory_cache.clear()
        if self.cache_file.exists():
            self.cache_file.unlink()
        self.is_dirty = False
        self.save_counter = 0
    
    def size(self) -> int:
        """Get current cache size."""
        return self.memory_cache.size()
    
    def keys(self):
        """Get cache keys."""
        return self.memory_cache.keys()
    
    def values(self):
        """Get cache values."""
        return self.memory_cache.values()
    
    def items(self):
        """Get cache items."""
        return self.memory_cache.items()


class CacheManager:
    """
    Centralized cache management for multiple named caches.
    
    Provides a registry of caches with unified management operations.
    """
    
    def __init__(self, cache_dir: Path):
        """
        Initialize cache manager.
        
        Parameters:
        -----------
        cache_dir : Path
            Directory for cache storage
        """
        self.cache_dir = cache_dir
        self.caches: Dict[str, PersistentCache] = {}
        
        # Ensure cache directory exists
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def register_cache(
        self,
        name: str,
        cache_file: Optional[Path] = None,
        **kwargs
    ) -> PersistentCache:
        """
        Register a new cache.
        
        Parameters:
        -----------
        name : str
            Cache name
        cache_file : Path, optional
            Path to cache file (defaults to cache_dir/name.json)
        **kwargs : dict
            Additional arguments for PersistentCache
        
        Returns:
        --------
        PersistentCache
            The registered cache instance
        """
        if cache_file is None:
            cache_file = self.cache_dir / f"{name}.json"
        
        cache = PersistentCache(cache_file, **kwargs)
        self.caches[name] = cache
        return cache
    
    def get_cache(self, name: str) -> Optional[PersistentCache]:
        """Get a cache by name."""
        return self.caches.get(name)
    
    def save_all(self) -> None:
        """Save all registered caches."""
        for name, cache in self.caches.items():
            cache.save()
            logger.debug(f"Saved cache: {name}")
    
    def clear_all(self) -> None:
        """Clear all registered caches."""
        for name, cache in self.caches.items():
            cache.clear()
            logger.debug(f"Cleared cache: {name}")
    
    def stats(self) -> Dict[str, int]:
        """
        Get statistics for all caches.
        
        Returns:
        --------
        dict
            Mapping of cache name to size
        """
        return {name: cache.size() for name, cache in self.caches.items()}
