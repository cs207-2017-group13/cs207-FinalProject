from setuptools import setup


setup(
    name="chemical_kinetics",
    version="0.0.1",
    author=["Divyam Misra", "Shiyun Qiu", "Victor Zhao"],
    description=("Chemical kinetics library"),
    license="BSD",
    package_dir={
        '': "src"
    },
    package_data={
        '': ["data"]
    }
)
