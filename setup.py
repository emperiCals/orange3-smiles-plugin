from setuptools import setup, find_packages

setup(
    name="Orange3-SMILES-Plugin",
    version="1.0.2",
    description="A plugin to extract molecular descriptors and visualize SMILES.",
    packages=find_packages(),
    install_requires=["Orange3", "rdkit", "numpy", "AnyQt"],
    entry_points={
        "orange.widgets": (
            "Data = orangecontrib.smiles_features.data",
            "Visualize = orangecontrib.smiles_features.visualize",
        ),
    },
)