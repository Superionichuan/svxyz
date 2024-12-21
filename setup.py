from setuptools import setup, find_packages

setup(
    name="svxyz",
    version="2.1.6",
    description="A suite of scripts for XYZ data processing (txyz, xyz2pos, dxyz, pxyz)",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Shichuan Sun",
    author_email="einstchuan@gmail.com",
    url="https://github.com/Superionichuan/svxyz",
    packages=find_packages(),
    install_requires=[
        "ase",
        "numpy",
        "matplotlib",
        "scipy",
        "mplcursors"
    ],
    entry_points={
        "console_scripts": [
            "txyz=svxyz.txyz:main",
            "xyz2pos=svxyz.xyz2pos:main",
            "dxyz=svxyz.dxyz:main",
            "pxyz=svxyz.pxyz:main",
            "analpos=svxyz.analpos:main",
            "asefmt=svxyz.asefmt:main"
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)

