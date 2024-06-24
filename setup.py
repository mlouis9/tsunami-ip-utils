from setuptools import setup, find_packages

setup(
    name="your_package_name",
    version="0.1.0",
    author="Matthew Louis",
    author_email="matthewlouis31@gmail.com",
    description="A tool for visualizing similarity data from the SCALE code: TSUNAMI-IP",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="http://github.com/mlouis9/tsunami_ip_utils",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        'plotly~=5.22.0',
        'matplotlib',
        'pandas',
        'uncertainties',
        'scipy',
        'flask',
        'dash'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)