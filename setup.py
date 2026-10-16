from setuptools import find_namespace_packages, setup


setup(
    name="dn-handler",
    version="0.1.0",
    description="Utilities for handling de-novo sequence data.",
    packages=find_namespace_packages(include=["dn_handler", "dn_handler.*"]),
    python_requires=">=3.8",
)
