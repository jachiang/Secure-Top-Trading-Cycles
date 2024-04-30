<!-- Install OpenFHE -->
<!-- Run: cmake .. -->
<!-- Run: make all -->

# Secure Top Trading Cycle Implementation

Implementation of the top trading cycle algorithm over packed ciphertexts.

- Install [OpenFHE](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html) and required dependencies.
- Run `cmake CMakeLists.txt` in repository to generate build files.
- Run `make all` in repository to compile `secure_cycle_finding.cpp`.
- Run `./secure_cycle_finding` to execute compiled benchmark binary.
- Set `numParties` to 5, 10, 15, 20, 25 in L65 of `secure_cycle_finding.cpp` to benchmark different number of parties.
