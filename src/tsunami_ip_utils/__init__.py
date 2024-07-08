import yaml
from pathlib import Path

current_dir = Path(__file__).resolve().parent
def load_config():
    with open(current_dir / "config.yaml", 'r') as stream:
        return yaml.safe_load(stream)

# Load the configuration when the package is imported
config = load_config()