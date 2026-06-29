# Contributor Guide

## Building the Docs

To build and preview the documentation locally, create a new Conda environment from the `docs-environment.yml` and activate it.
Then, `cd` into the `doc/` directory and run `make html`.

```sh
conda env create -f docs-environment.yml -n docs
conda activate docs
cd doc
make html
```

You can view the documentation in your browser by opening `docs/build/html/index.html`.

## Creating Static Pages

Simple pages may be written in Markdown  or [reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) and stored in the `doc/source/` directory.
Markdown support is enabled via the MyST parser extension, refer to the [MyST syntax guide](https://myst-parser.readthedocs.io/en/latest/syntax/typography.html)
for more information related to composition of Markdown documents.
Add the file to a `toc` block in the `doc/source/index.rst` in order for it be included in the main documentation.


1. Create your document containing useful information under `doc/source/useful_info.md`.

```md
<!-- doc/source/useful_info.md -->

# Useful Information
```

2. Add the file to a `toc` block in the `doc/source/index.rst`.

```rst
.. toctree::

   useful_info
```

:::{tip}
If you want to create a hierarchical structure of pages, add the subdirectory to the `doc/source/` directory and add the file to a `toc` block in the parent file.
:::

## Documenting Python Code

### Where to add documentation

Add or edit docstrings directly in the Python source files following the NumPy docstring style.
For generated API pages, autosummary will produce stub pages from the docstrings and signatures.
Any new top-level `orbit` submodule must be registered by adding it to `docs/source/modules.rst`.

### Docstring format

Use the [NumPy docstring convention](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard).
Keep short summary on the first line.
Use Parameters with a table-like list; include types and optional default values.

```python
def foo(x, y=1):
    """
    Short one-line summary.

    Longer description (optional), explaining behavior, side effects, and important details.

    .. math::

        z = x \cdot y

    Parameters
    ----------
    x : int
        Description of x.
    y : int, optional
        Description of y. Defaults to 1.

    Returns
    -------
    float
        Description of the return value.

    Raises
    ------
    ValueError
        If x is negative.

    """
    ...
```

## Documenting C++ Code

### Where to add documentation

Add Doxygen comments directly in the C++ header or source files above the declarations.

Example documenting a function:

```cpp
/**
 * @brief Compute a simple value.
 *
 * The relation is given by the LaTeX expression: \f\$ z = x \cdot y \f\$.
 * Detailed description if needed.
 *
 * @param x First input.
 * @param y Second input.
 * @return Computed value.
 * @throws std::invalid\_argument if x <= 0.
 *
 * Example:
 * @code
 * double r = foo(x, y);
 * @endcode
 */
double foo(double x, double y);
```

Example documenting a class:

```cpp
/**
 * @class Bar
 * @brief Represents a simple concept.
 *
 * More detailed class description.
 *
 * Mathematical note: \f\$E = mc^2\f\$
 */
class Bar {
public:
    /** Construct with a value. */
    Bar(double v);
    /** Return the value. */
    double value() const;
    // ...
};
```

## For Maintainers

The docs are deployed to the [main GitHub Pages repository for the PyORBIT Collaboration](https://github.com/pyorbit-collaboration/pyorbit-collaboration.github.io)
via an Action ([peaceiris/actions-gh-pages](https://github.com/peaceiris/actions-gh-pages)) defined in [PyORBIT-Collaboration/PyORBIT3/.github/workflow/docs.yml](https://github.com/PyORBIT-Collaboration/PyORBIT3/blob/main/.github/workflows/docs.yml).
The action requires a [Deploy Key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/managing-deploy-keys) to be added to the PyORBIT-Collaboration/PyORBIT-Collaboration.github.io repository
as well as a [Secret](https://docs.github.com/en/actions/security-guides/automatic-token-authentication) to be created in the PyORBIT-Collaboration/PyORBIT3 repository.
Refer to the [instructions for deploying to an external repository](https://github.com/peaceiris/actions-gh-pages?tab=readme-ov-file#%EF%B8%8F-deploy-to-external-repository-external_repository) and [Create SSH Deploy Key](https://github.com/peaceiris/actions-gh-pages?tab=readme-ov-file#tips-and-faq) for more details.
