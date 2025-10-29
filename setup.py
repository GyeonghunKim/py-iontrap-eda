import io
from setuptools import find_packages, setup


# Read in the README for the long description on PyPI
def long_description():
    with io.open("README.md", "r", encoding="utf-8") as f:
        readme = f.read()
    return readme


setup(
    name="py_iontrap_eda",
    version="0.3",
    description="Ion trap chip EDA tool",
    long_description=long_description(),
    url="https://github.com/gyeonghun-kim/py_iontrap_eda",
    author="Gyeonghun Kim",
    author_email="gyeonghun.kim@duke.edu",
    license="MIT",
    packages=find_packages(),
    classifiers=["Programming Language :: Python :: 3"],
    zip_safe=False,
)
