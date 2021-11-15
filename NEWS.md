# morgancpp 0.4.0

* New `MorganMap` data structure for fast identity matching of fingerprints
* Support for reading RDkit run-length encoded (RLE) fingerprints generated using
  `GetMorganFingerprintAsBitVect(inchi).ToBinary().hex()`.

# morgancpp 0.3.0

* New binary format which works for extremely large fingerprint vectors

# morgancpp 0.2.0

* New binary format for fingerprints on disk incompatible with previous versions
* Fast loading and reading using Zstandard compression
