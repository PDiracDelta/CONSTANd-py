# INSTALL INSTRUCTIONS:
After cloning this repository and entering the top directory, use the following command to install in-place, in editable mode:
`python3 -m pip install -e <path_to_repository>`
You can update it by simply pulling from the repository again.

# USAGE
To then use it in Python3, simply import it and use it on some matrix-like data:
```
from constand.constand import constand
# generate some data
import numpy as np
m = np.random.rand(100, 6)
# normalize it
result = constand(m)
normalized_m = result['normalizedData']

# available metadata
result.keys()
```

> dict_keys(['normalizedData', 0, 'convergenceTrail', 1, 'R', 2, 'S', 3])

The numerical entries are for backwards compatibility (the result used to be a list).
