from importlib import import_module
from pathlib import Path

commands = {}

for f in sorted(Path(__file__).parent.glob('*.py')):
    if '__' in f.stem:
        continue
    commands[f.stem.replace('_', '-')] = import_module(f'.{f.stem}', __package__)
