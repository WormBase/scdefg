import setuptools

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    setuptools.setup(
        name="scdefg",
        version="0.2.0",
        author="Eduardo da Veiga Beltrame",
        author_email="munfred@brandeis.edu",
        description="A single page Flask app with GUI for performing differential expression on with scvi-tools.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/munfred/scdefg",
        license='LICENSE.txt',
        packages=setuptools.find_packages(),
        install_requires=requirements,
        classifiers=[
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "License :: OSI Approved :: The Unlicense (Unlicense)",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
        python_requires='>=3.6',
    )

