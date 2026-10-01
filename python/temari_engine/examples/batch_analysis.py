"""Small end-to-end examples for the published NumPy batch readers.

Run after installing this package with Python 3.10+:
    python python/temari_engine/examples/batch_analysis.py
"""

import numpy as np

from temari_engine import data as td


def main():
    s_values = np.linspace(0.0, 6.0, 1001)  # Angstrom^-1

    f = td.open("F")
    f_batch = f.channel(79, "L3").f_many(200.0, s_values)
    print("Au L3 F:", f_batch.values.shape, sorted(set(f_batch.regions)))

    neutral_ds = td.open("factors")
    neutral_elements = [neutral_ds.element(z, tag) for z, tag in neutral_ds.species()]
    neutral = neutral_ds.evaluate_many(neutral_elements, "fx", s_values)
    print("neutral f_x:", neutral.values.shape, neutral.species[25])

    ion_ds = td.open("factors-ion")
    ion_elements = [ion_ds.element(z, tag) for z, tag in ion_ds.species()]
    ions = ion_ds.evaluate_many(ion_elements, "fe_regular", s_values)
    print("ion f_e regular:", ions.values.shape, ions.species[0], ions.source.citation())


if __name__ == "__main__":
    main()
