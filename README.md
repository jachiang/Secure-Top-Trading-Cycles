# Secure Top Trading Cycles Implementation

Implementation of the top trading cycles algorithm over packed BGW ciphertexts.

- Install [OpenFHE](https://github.com/openfheorg/openfhe-hexl) with HEXL back-end (requires Intel AVX-512 support).
- Run `cmake CMakeLists.txt` in repository to generate build files.
- Run `make all` in repository to compile `secure_cycle_finding.cpp`.
- Run `./secure_cycle_finding` to execute compiled benchmark binary.
- Set `numParties` to 5, 10, 15, 20, 25 in L65 of `secure_cycle_finding.cpp` to benchmark different number of parties.
