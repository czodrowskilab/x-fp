# Contributing to X-FP

## Introduction
Thanks for your interest in contributing to X-FP! We welcome contributions from everyone. Please read this document before you start contributing.


## What can I contribute?
We welcome contributions of all kinds, including but not limited to:

1. Bug reports
2. Bug fixes
3. New features
4. Documentation improvements
5. Code quality improvements
6. Test improvements
7. Examples

We will expand on these topics below.

### Bug reports
Something not working as expected? Please open an issue on GitHub. Kindly include the following information in your bug report:

- A short description of the bug.
- A minimal code snippet to reproduce the bug.
- The expected behavior.
- The actual behavior.
- The version of X-FP you are using.
- The versions of the dependencies you are using.

### Bug fixes
Squashed a bug? So, what's next? Please open a pull request on GitHub. 

### New features
X-FP is an organic project and we are always looking for new features. 
However, some of the feature addition we would like to see are:

- Implementation of additional molecular fingerprints.
- Adding new feature importance methods (see tutorial on that [here](https://github.com/czodrowskilab/x-fp/blob/contributing/adding_feature_importance_method.md)).
- Adding new X-FP output implementations, for instance, Jupyter notebook output, machine-readable outputs, etc.
- Single molecule analysis.

Interested in seeing a new feature in X-FP but can not do it yourself at the moment? No worries! Just open an issue on GitHub and/or contact us.

### Documentation improvements
X-FP uses [Sphinx](https://www.sphinx-doc.org/en/master/) for documentation. The docstrings are written in [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) style.

### Code quality improvements
Think the code can be improved? Please open a pull request on GitHub. We welcome contributions to improve the code quality.

### Test improvements
X-FP uses [pytest](https://docs.pytest.org/en) for testing. The tests are located in the `tests` folder. So far, we only test the core functionality of X-FP (`core_xfp.py`). We welcome contributions to expand the test coverage.

### Examples
Have a cool example of how you used X-FP and would like to share it with the community? Please do so!


## How can I contribute?

If you want to add a new feature, please open a pull request on GitHub. Kindly include the following information in your pull request:

- A short description of the new feature.
- The reason why you think this feature should be added to X-FP with an example use case.
- A minimal code snippet to demonstrate the new feature.
- A short description of the test you have added to test the new feature.
- Documentation updates, if applicable.
- Run `pytest` to make sure all tests pass. If you have added new tests, please make sure they pass as well.

X-FP uses [black](https://github.com/psf/black) for code formatting. Please make sure your code is formatted using `black` before opening a pull request. You must install `black` in your local environment to do this. Then, simply run `black .` in the root directory of X-FP to format the code. The `black` configurations are already present in `pyproject.toml` file. 

Adapt above steps as necessary for other types of contributions.