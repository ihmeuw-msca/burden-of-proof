import importlib
from pathlib import Path
from setuptools import setup, find_packages


if __name__ == "__main__":
    base_dir = Path(__file__).parent
    spec = importlib.util.spec_from_file_location(
        "__about__",
        base_dir / "src" / "espipeline" / "__about__.py"
    )
    about = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(about)

    with (base_dir/"README.md").open() as f:
        long_description = f.read()

    install_requirements = [
        "numpy",
        "scipy",
        "pandas",
    ]

    test_requirements = [
        "pytest",
        "pytest-mock",
    ]

    doc_requirements = []

    setup(name=about.__title__,
          version=about.__version__,

          description=about.__summary__,
          long_description=long_description,
          license=about.__license__,
          url=about.__uri__,

          author=about.__author__,
          author_email=about.__email__,

          package_dir={"": "src"},
          packages=find_packages(where="src"),
          include_package_data=True,

          install_requires=install_requirements,
          tests_require=test_requirements,
          extras_require={
              "docs": doc_requirements,
              "test": test_requirements,
              "dev": doc_requirements + test_requirements
          },
          zip_safe=False,)
