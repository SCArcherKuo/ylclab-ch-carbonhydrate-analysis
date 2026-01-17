"""
Retry Utilities

This module provides retry decorators and utilities for handling transient failures
in API calls with exponential backoff.
"""

import time
import functools
from typing import Callable, TypeVar, Any, Tuple, Type
from loguru import logger

T = TypeVar('T')


def retry_with_backoff(
    max_retries: int = 3,
    base_delay: float = 2.0,
    backoff_factor: float = 2.0,
    exceptions: Tuple[Type[Exception], ...] = (Exception,),
    on_retry: Callable[[Exception, int], None] = None
):
    """
    Decorator that retries a function with exponential backoff.
    
    Parameters:
    -----------
    max_retries : int
        Maximum number of retry attempts
    base_delay : float
        Initial delay between retries in seconds
    backoff_factor : float
        Multiplier for delay after each retry
    exceptions : tuple of Exception types
        Exception types to catch and retry
    on_retry : callable, optional
        Callback function called on each retry with (exception, attempt_number)
    
    Returns:
    --------
    decorator
        Function decorator that adds retry logic
    
    Examples:
    ---------
    >>> @retry_with_backoff(max_retries=3, base_delay=1.0)
    ... def unstable_api_call():
    ...     return requests.get(url)
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> T:
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    if attempt == max_retries - 1:
                        # Last attempt failed, re-raise
                        raise
                    
                    # Calculate delay with exponential backoff
                    delay = base_delay * (backoff_factor ** attempt)
                    
                    # Call retry callback if provided
                    if on_retry:
                        on_retry(e, attempt + 1)
                    else:
                        logger.debug(
                            f"Retry {attempt + 1}/{max_retries} after {delay:.1f}s "
                            f"for {func.__name__}: {str(e)}"
                        )
                    
                    time.sleep(delay)
            
            # Should never reach here
            raise RuntimeError(f"Retry logic failed for {func.__name__}")
        
        return wrapper
    return decorator


class RetryManager:
    """
    Manages retry logic with configurable parameters.
    
    This class provides a more flexible approach to retries than the decorator,
    allowing for runtime configuration and state tracking.
    """
    
    def __init__(
        self,
        max_retries: int = 3,
        base_delay: float = 2.0,
        backoff_factor: float = 2.0
    ):
        """
        Initialize retry manager.
        
        Parameters:
        -----------
        max_retries : int
            Maximum number of retry attempts
        base_delay : float
            Initial delay between retries in seconds
        backoff_factor : float
            Multiplier for delay after each retry
        """
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.backoff_factor = backoff_factor
        self.retry_count = 0
    
    def execute_with_retry(
        self,
        func: Callable[..., T],
        *args,
        exceptions: Tuple[Type[Exception], ...] = (Exception,),
        on_retry: Callable[[Exception, int], None] = None,
        **kwargs
    ) -> T:
        """
        Execute a function with retry logic.
        
        Parameters:
        -----------
        func : callable
            Function to execute
        *args : tuple
            Positional arguments for func
        exceptions : tuple of Exception types
            Exception types to catch and retry
        on_retry : callable, optional
            Callback function called on each retry
        **kwargs : dict
            Keyword arguments for func
        
        Returns:
        --------
        any
            Result of func execution
        
        Raises:
        -------
        Exception
            Re-raises the last exception if all retries fail
        """
        for attempt in range(self.max_retries):
            try:
                result = func(*args, **kwargs)
                self.retry_count = 0  # Reset on success
                return result
            except exceptions as e:
                if attempt == self.max_retries - 1:
                    raise
                
                delay = self.base_delay * (self.backoff_factor ** attempt)
                
                if on_retry:
                    on_retry(e, attempt + 1)
                else:
                    logger.debug(
                        f"Retry {attempt + 1}/{self.max_retries} after {delay:.1f}s: {str(e)}"
                    )
                
                time.sleep(delay)
                self.retry_count += 1
        
        raise RuntimeError("Retry logic failed unexpectedly")
    
    def reset(self):
        """Reset retry counter."""
        self.retry_count = 0
