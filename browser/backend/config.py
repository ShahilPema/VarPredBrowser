"""
Configuration loader for VarPredBrowser

Loads configuration from config/paths.yaml with support for environment variable overrides.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional


def get_project_root() -> Path:
    """Get the project root directory (where config/ lives)."""
    # Start from this file's location and go up to find project root
    current = Path(__file__).resolve().parent  # browser/backend/

    # Go up until we find config/paths.yaml
    for _ in range(5):  # Max 5 levels up
        if (current / 'config' / 'paths.yaml').exists():
            return current
        current = current.parent

    # Fall back to working directory
    return Path.cwd()


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file.

    Args:
        config_path: Optional path to config file. If not provided,
                     looks for config/paths.yaml in project root.

    Returns:
        Configuration dictionary with all settings.
    """
    if config_path is None:
        project_root = get_project_root()
        config_path = project_root / 'config' / 'paths.yaml'
    else:
        config_path = Path(config_path)

    if not config_path.exists():
        print(f"Warning: Config file not found at {config_path}, using defaults")
        return get_default_config()

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Ensure browser section exists with defaults
    if 'browser' not in config:
        config['browser'] = {}

    # Apply browser defaults
    browser_defaults = {
        'data_dir': 'data',
        'bigwig_dir': '../GenomeBrowser2/test_data/BW',
        'structures_dir': 'data/structures',
        'host': '0.0.0.0',
        'port': 8000,
    }

    for key, default_value in browser_defaults.items():
        if key not in config['browser']:
            config['browser'][key] = default_value

    # Apply environment variable overrides
    apply_env_overrides(config)

    return config


def apply_env_overrides(config: Dict[str, Any]) -> None:
    """Apply environment variable overrides to config."""

    # Browser data directory
    if os.getenv('VARPRED_DATA_DIR'):
        config['browser']['data_dir'] = os.getenv('VARPRED_DATA_DIR')

    # BigWig directory
    if os.getenv('VARPRED_BIGWIG_DIR'):
        config['browser']['bigwig_dir'] = os.getenv('VARPRED_BIGWIG_DIR')

    # Server settings
    if os.getenv('VARPRED_HOST'):
        config['browser']['host'] = os.getenv('VARPRED_HOST')

    if os.getenv('VARPRED_PORT'):
        config['browser']['port'] = int(os.getenv('VARPRED_PORT'))


def get_default_config() -> Dict[str, Any]:
    """Return default configuration."""
    return {
        'browser': {
            'data_dir': 'data',
            'bigwig_dir': '../GenomeBrowser2/test_data/BW',
            'structures_dir': 'data/structures',
            'host': '0.0.0.0',
            'port': 8000,
        },
        'data_sources': {},
        'output': {
            'curve_data': 'output/curve_data',
            'figures': 'output/figures',
        }
    }


def resolve_path(path: str, base_dir: Optional[Path] = None) -> Path:
    """
    Resolve a path, making it absolute if relative.

    Args:
        path: Path string (can be absolute or relative)
        base_dir: Base directory for relative paths (defaults to project root)

    Returns:
        Resolved absolute Path
    """
    p = Path(path)
    if p.is_absolute():
        return p

    if base_dir is None:
        base_dir = get_project_root()

    return (base_dir / p).resolve()


# Global config instance (loaded on first access)
_config: Optional[Dict[str, Any]] = None


def get_config() -> Dict[str, Any]:
    """Get the global configuration, loading it if necessary."""
    global _config
    if _config is None:
        _config = load_config()
    return _config


def get_data_dir() -> Path:
    """Get the browser data directory."""
    config = get_config()
    return resolve_path(config['browser']['data_dir'])


def get_bigwig_dir() -> Path:
    """Get the BigWig files directory."""
    config = get_config()
    return resolve_path(config['browser']['bigwig_dir'])


def get_structures_dir() -> Path:
    """Get the structures directory."""
    config = get_config()
    return resolve_path(config['browser']['structures_dir'])
