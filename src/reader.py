import numpy as np
import pandas as pd
import zstandard
import json

def open_posteriors(filename):
    meta = json.load(open(filename + ".meta"))
    raw_floats = np.frombuffer(
        zstandard.ZstdDecompressor().decompress(
            open(filename, "rb").read(),
            max_output_size = meta["sequence_length"] * meta["num_pairs"] * 2 * 4,
        ),
        dtype=np.float32
    )

    n_chunks = int(np.ceil(meta["num_pairs"] / meta["chunk_size"]))

    raw_floats = raw_floats.reshape((
        n_chunks,
        2,
        meta["sequence_length"],
        meta["chunk_size"],
    ))
    
    raw_floats = np.transpose(raw_floats, [1, 2, 0, 3]).reshape(
            (2, meta["sequence_length"], -1)
        )[:, :, :meta["num_pairs"]]
    
    column_names = [f"{i}_{j}" for i,j in meta["pairs"]]
    alphas = pd.DataFrame(raw_floats[0, :, :], index=meta["output_positions"], columns=column_names)
    betas = pd.DataFrame(raw_floats[1, :, :], index=meta["output_positions"], columns=column_names)

    return alphas, betas, meta
    





