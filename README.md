# Meteor: Improved Secure 3-Party Neural Network Inference with Reducing Online Communication Costs

A privacy-preserving framework for efficient 3-party protocols tailored for neural networks, the paper is avaiable [eprint/2023/100](https://eprint.iacr.org/2023/100). This work builds on [Falcon](https://snwagh.github.io).


### Table of Contents

- [Warning](#warning)
- [Requirements](#requirements)
- [Source Code](#source-code)
    - [Repository Structure](#repository-structure)
    - [Building the code](#building)
    - [Running the code](#running)


### Warning
---
This codebase is released solely as a reference for other developers, as a proof-of-concept, and for benchmarking purposes. In particular, it has not had any security review, has a number of implementational TODOs, has a number of known bugs (especially in the malicious implementation), and thus, should be used at your own risk. You can contribute to this project by creating pull requests and submitting fixes and implementations. The code has not run end-to-end inference and we provide online computation at present, the setup phase will be added in future.


### Requirements
---
* The code should work on most Linux distributions (It has been developed and tested with [Ubuntu](http://www.ubuntu.com/) 18.04).

* **Required packages for Meteor:**
  * [`g++`](https://packages.debian.org/testing/g++)
  * [`make`](https://packages.debian.org/testing/make)
  * [`libssl-dev`](https://packages.debian.org/testing/libssl-dev)

  Install these packages with your favorite package manager, e.g, `sudo apt-get install <package-name>`.

### Source Code
---

#### Repository Structure

* `files/`    - Shared keys, IP addresses and data files.
* `files/preload`    - Contains data for pretrained network from SecureML. The other networks can be generated using `scripts` and functions in `secondary.cpp`
* `lib_eigen/`    - [Eigen library](http://eigen.tuxfamily.org/) for faster matrix multiplication.
* `src/`    - Source code.
* `util/` - Dependencies for AES randomness.
* `scripts/` - Contains python code to generate trained models for accuracy testing over a batch.
* The `god` script makes remote runs simpler (as well as the `makefile`)

#### Building the code

To build Meteor, run the following commands:

```
cd Meteor
make all -j$(nproc)
```

#### Running the code

To run the code, simply choose one of the following options: 

* `make`: Prints all the possible makefile options.
* `make terminal`: Runs the 3PC code on localhost with output from $P_0$ printed to standard output.
* `make file`: : Runs the 3PC code on localhost with output from $P_0$ printed to a file (in `output/3PC.txt`)
* `make valg`: Useful for debugging the code for set faults. Note that the -03 optimization flag needs to be suppressed (toggle lines 42, 43 in `makefile`)
* `make command`: Enables running a specific network, dataset, adversarial model, and run type (localhost/LAN/WAN) specified through the `makefile`. This takes precedence over choices in the `src/main.cpp` file.
* To run the code over tmux over multiple terminals, `make zero`, `make one`, and `make two` come in handy.
* Finally, the `makefile` (line 4-15) contains the descriptions of the arguments accepted by the executable.

---
For questions, please create git issues; for eventual replies, you can also reach out to [19950512dy@gmail.com](19950512dy@gmail.com)
