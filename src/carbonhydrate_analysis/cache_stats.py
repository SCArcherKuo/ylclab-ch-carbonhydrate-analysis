"""
Cache Statistics Utility

Display statistics about the ChEBI cache.
"""

from pathlib import Path
import json


def show_cache_stats():
    """Display statistics about the ChEBI cache."""
    project_root = Path(__file__).parent.parent.parent
    cache_file = project_root / 'data' / 'cache' / 'chebi_children_cache.json'
    
    if not cache_file.exists():
        print("No cache file found.")
        print(f"Expected location: {cache_file}")
        return
    
    # Load cache
    with open(cache_file, 'r', encoding='utf-8') as f:
        cache_data = json.load(f)
    
    # Calculate statistics
    total_entries = len(cache_data)
    total_children = sum(len(children) for children in cache_data.values())
    file_size = cache_file.stat().st_size
    
    # Find entries with most children
    entries_by_children = sorted(
        [(chebi_id, len(children)) for chebi_id, children in cache_data.items()],
        key=lambda x: x[1],
        reverse=True
    )
    
    print("=" * 60)
    print("ChEBI Cache Statistics")
    print("=" * 60)
    print(f"\nCache file: {cache_file}")
    print(f"File size: {file_size:,} bytes ({file_size / 1024 / 1024:.2f} MB)")
    print(f"\nTotal ChEBI IDs cached: {total_entries}")
    print(f"Total children relationships: {total_children}")
    print(f"Average children per ID: {total_children / total_entries:.1f}")
    
    print(f"\n{'ChEBI ID':<15} {'Children Count':<15}")
    print("-" * 30)
    for chebi_id, count in entries_by_children[:10]:
        print(f"{chebi_id:<15} {count:<15}")
    
    if len(entries_by_children) > 10:
        print(f"... and {len(entries_by_children) - 10} more")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    show_cache_stats()
