from setuptools import setup, find_packages

def read_requirements(file):
    with open(file) as f:
        return f.read().splitlines()

def read_file(file):
   with open(file) as f:
        return f.read()
    
long_description = read_file("README.rst")
version = read_file("VERSION")
requirements = read_requirements("requirements.txt")

setup(
    name="ionic_strength_calc_of_tris-EDTA",
    version="0.1.0",
    packages=["ionic_strength_calc_of_tris-EDTA",],
    install_requires=[],
    license="MIT",
    url="https://github.com/oskarestrom/nd2handling/ionic_strength_calc_of_tris-EDTA",
    author="Oskar Ström",
    author_email="oskarestrom@protonmail.com",
    description="Calculates the Ionic Strength of a Tris-EDTA solution",
    long_description=long_description,
    long_description_content_type="text/x-rst", # If this causes a warning, upgrade your setuptools package
    packages = find_packages(exclude=["test"]),  # Don't include test directory in binary distribution
    install_requires = requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]  # Update these accordingly
)

