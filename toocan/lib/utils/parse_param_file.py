import re
from collections import defaultdict

def parse_param_file(filepath):
    params = defaultdict(list)

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip blank lines and comments
            if not line or line.startswith(';') or '=' not in line:
                continue

            # Parse key and value
            key, value = map(str.strip, line.split('=', 1))

            # Try converting value to int, then float, then leave as string
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    value = value.strip('"').strip("'")  # remove quotes if any

            params[key].append(value)

    # Simplify lists that only have one item
    clean_params = {k: v[0] if len(v) == 1 else v for k, v in params.items()}
    return clean_params
