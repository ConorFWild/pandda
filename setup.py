import setuptools
from setuptools.command.install import install


with open("README.md", "r") as fh:
    long_description = fh.read()

def requirements():
    # The dependencies are the same as the contents of requirements.txt
    with open('requirements.txt') as f:
        return [line.strip() for line in f if line.strip()]

setuptools.setup(
    name="pandda_2",
    version="0.0.3",
    author="Conor Francis Wild",
    author_email="conor.wild@sky.com",
    cmdclass={"install": CustomInstallCommand},
    description="A package for handling many crystalographic datasets simultainiously",
    long_description="",
    long_description_content_type="text/markdown",
    url="https://github.com/ConorFWild/pandda.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=["program/pandda.analyse",
             "program/pandda.inspect",
             ],
    python_requires='>=2.7',
    install_requires=requirements(),
)
