from setuptools import setup, find_packages


setup(
    name="chemical-kinetics",
    version="0.0.1",
    author=["Divyam Misra", "Shiyun Qiu", "Victor Zhao"],
    description=("Chemical kinetics library"),
    license="BSD",
    packages=find_packages('src'),
    package_dir={"": "src"},
    package_data={'chemkin': ["*.sqlite"]},
    scripts=['src/bin/process_reaction_system'],
    install_requires=["scipy"]
)
